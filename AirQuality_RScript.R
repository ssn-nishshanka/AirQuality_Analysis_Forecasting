
library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)
library(parallel)
library(ggplot2)

data <- read_csv(file.choose())

glimpse(data)
summary(data)

str(data)

sum(is.na(data$co))

median_co <- median(data$co, na.rm = TRUE)
data <- data %>%
  mutate(co = ifelse(is.na(co), median_co, co))  

boxplot(data$co, main = "Boxplot of CO Levels", ylab = "CO")

Q1 <- quantile(data$co, 0.25, na.rm = TRUE)
Q3 <- quantile(data$co, 0.75, na.rm = TRUE)
IQR_value <- IQR(data$co, na.rm = TRUE)

lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value 

data <- data %>%
  mutate(co = ifelse(co < lower_bound | co > upper_bound, median_co, co))

monthly_data <- data %>%
  mutate(month = floor_date(tanggal, "month")) %>%
  group_by(month) %>%
  summarise(monthly_co = mean(co, na.rm = TRUE))

head(monthly_data)

ggplot(monthly_data, aes(x = month, y = monthly_co)) +
  geom_line(color = "blue") +
  labs(title = "Monthly CO Levels", x = "Month", y = "Monthly CO Level (µg/m³)") +
  theme_minimal()

monthly_ts <- ts(monthly_data$monthly_co, frequency = 12, start = c(2010, 1))

decomposed_ts_stl <- stl(monthly_ts, s.window = "periodic")

plot(decomposed_ts_stl)

trend <- decomposed_ts_stl$time.series[, "trend"]
seasonal <- decomposed_ts_stl$time.series[, "seasonal"]
residuals <- decomposed_ts_stl$time.series[, "remainder"]

plot(trend, main = "Trend")
plot(seasonal, main = "Seasonal")
plot(residuals, main = "Residuals")

acf(monthly_ts, main = "Autocorrelation of Monthly CO Levels")

pacf(monthly_ts, main = "Partial Autocorrelation of Monthly CO Levels")

mean_co <- mean(monthly_data$monthly_co, na.rm = TRUE)
variance_co <- var(monthly_data$monthly_co, na.rm = TRUE)

cat("Mean of CO Levels:", mean_co, "\n")
cat("Variance of CO Levels:", variance_co, "\n")

adf_result <- adf.test(monthly_data$monthly_co, alternative = "stationary")

print(adf_result)

monthly_diff <- diff(monthly_ts, differences = 1)

par(mfrow = c(1, 2))
acf(monthly_diff, main = "ACF of Differenced Monthly CO Data")
pacf(monthly_diff, main = "PACF of Differenced Monthly CO Data")
par(mfrow = c(1, 1))

adf_result_diff <- adf.test(monthly_diff, alternative = "stationary")
print(adf_result_diff)

train_size <- round(length(monthly_diff) * 0.8)
train_data <- monthly_diff[1:train_size]

test_data <- monthly_diff[(train_size + 1):length(monthly_diff)]

train_start_year <- start(monthly_diff)[1] + (start(monthly_diff)[2] - 1) / 12
train_data <- ts(train_data, frequency = 12, start = c(train_start_year))

fit_arima <- arima(train_data, order = c(3, 1, 1)) 

summary(fit_arima)

checkresiduals(fit_arima)

monthly_diff_seasonal <- diff(monthly_diff, lag = 12)

train_size_seasonal <- round(length(monthly_diff_seasonal) * 0.8)
train_data_seasonal <- monthly_diff_seasonal[1:train_size_seasonal]

test_data_seasonal <- monthly_diff_seasonal[(train_size_seasonal + 1):length(monthly_diff_seasonal)]

train_start_year_seasonal <- start(monthly_diff_seasonal)[1] + (start(monthly_diff_seasonal)[2] - 1) / 12
train_data_seasonal <- ts(train_data_seasonal, frequency = 12, start = c(train_start_year_seasonal))

par(mfrow = c(1, 2))
acf(train_data_seasonal, main = "Seasonal ACF of Differenced Data")
pacf(train_data_seasonal, main = "Seasonal PACF of Differenced Data")
par(mfrow = c(1, 1))

adf_result_seasonal <- adf.test(monthly_diff_seasonal, alternative = "stationary")
print(adf_result_seasonal)

fit_sarima <- arima(train_data_seasonal, order = c(3, 1, 1), seasonal = list(order = c(2, 1, 1), period = 12))

summary(fit_sarima)

checkresiduals(fit_sarima)

fit_ets <- ets(train_data_seasonal)

summary(fit_ets)

checkresiduals(fit_ets)

arima_forecast <- forecast(fit_arima, h = length(test_data))
sarima_forecast <- forecast(fit_sarima, h = length(test_data))
ets_forecast <- forecast(fit_ets, h = length(test_data))

arima_accuracy <- accuracy(arima_forecast, test_data)
sarima_accuracy <- accuracy(sarima_forecast, test_data)
ets_accuracy <- accuracy(ets_forecast, test_data)

print(arima_accuracy)
print(sarima_accuracy)
print(ets_accuracy)

test_data_ts <- ts(test_data_seasonal, frequency = 12, start = start(train_data_seasonal) + c(0, length(train_data_seasonal) / 12))

sarima_forecast <- forecast(fit_sarima, h = length(test_data_seasonal))

plot(sarima_forecast, main = "SARIMA Forecast vs Actual Data", xlab = "Month", ylab = "CO Levels")
lines(test_data_ts, col = "red", lwd = 2) 

years_to_forecast <- 5
frequency <- 12
horizon <- years_to_forecast * frequency

sarima_5yr_forecast <- forecast(fit_sarima, h = horizon)

plot(sarima_5yr_forecast, main = "5 Year SARIMA Forecast", xlab = "Month", ylab = "CO Levels")
lines(test_data_ts, col = "red", lwd = 2)

accuracy_metrics <- accuracy(sarima_forecast, test_data_seasonal)

print(accuracy_metrics)

mae <- accuracy_metrics["Test set", "MAE"]
rmse <- accuracy_metrics["Test set", "RMSE"]
mape <- accuracy_metrics["Test set", "MAPE"]

cat("Mean Absolute Error (MAE):", mae, "\n")
cat("Root Mean Squared Error (RMSE):", rmse, "\n")
cat("Mean Absolute Percentage Error (MAPE):", mape, "\n")

