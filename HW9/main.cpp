//
//  main.cpp
//  HW9
//
//  Created by Jake Hu on 12/3/21.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <chrono>

using namespace std;


double up_factor, uptick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, R, barrier_price;
int no_of_trials, no_of_divisions, no_of_sampling_instants;

unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator;

double get_uniform()
{
    // http://www.cplusplus.com/reference/random/exponential_distribution/
    uniform_real_distribution<double> distribution(0.0, 1.0);
    double number = distribution(generator);
    return (number);
}

double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a = fabs(z);
    double t = 1.0 / (1.0 + a * p);
    double b = c2 * exp((-z) * (z / 2.0));
    double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
    n = 1.0 - b * n;
    if (z < 0.0) n = 1.0 - n;
    return n;
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
    const double& K,       // strike (exercise) price,
    const double& r,       // interest rate
    const double& sigma,   // volatility
    const double& time)      // time to maturity
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
    double d2 = d1 - (sigma * time_sqrt);
    return S * N(d1) - K * exp(-r * time) * N(d2);
};

double closed_form_down_and_out_european_call_option()
{
    // I took this formula from Wilmott, Howison and Dewynne, "The Mathematics of Financial Derivatives"
    double K = (2 * risk_free_rate) / (volatility * volatility);
    double A = option_price_call_black_scholes(initial_stock_price, strike_price,
        risk_free_rate, volatility, expiration_time);
    double B = (barrier_price * barrier_price) / initial_stock_price;
    double C = pow(initial_stock_price / barrier_price, -(K - 1));
    double D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
    return (A - D * C);
}

double closed_form_down_and_in_european_put_option()
{
    // just making it easier by renaming the global variables locally
    double S = initial_stock_price;
    double r = risk_free_rate;
    double T = expiration_time;
    double sigma = volatility;
    double H = barrier_price;
    double X = strike_price;

    // Took these formulae from some online reference
    double lambda = (r + ((sigma * sigma) / 2)) / (sigma * sigma);
    double temp = 2 * lambda - 2.0;
    double x1 = (log(S / H) / (sigma * sqrt(T))) + (lambda * sigma * sqrt(T));
    double y = (log(H * H / (S * X)) / (sigma * sqrt(T))) + (lambda * sigma * sqrt(T));
    double y1 = (log(H / S) / (sigma * sqrt(T))) + (lambda * sigma * sqrt(T));
    return (-S * N(-x1) + X * exp(-r * T) * N(-x1 + sigma * sqrt(T)) +
        S * pow(H / S, 2 * lambda) * (N(y) - N(y1)) -
        X * exp(-r * T) * pow(H / S, temp) * (N(y - sigma * sqrt(T)) - N(y1 - sigma * sqrt(T))));
}

double option_price_put_black_scholes(const double& S,      // spot price
    const double& K,      // Strike (exercise) price,
    const double& r,      // interest rate
    const double& sigma,  // volatility
    const double& time)
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
    double d2 = d1 - (sigma * time_sqrt);
    return K * exp(-r * time) * N(-d2) - S * N(-d1);
};

double closed_form_down_and_out_european_put_option()
{
    double vanilla_put = option_price_put_black_scholes(initial_stock_price, strike_price,
        risk_free_rate, volatility, expiration_time);
    double put_down_in = closed_form_down_and_in_european_put_option();
    return (vanilla_put - put_down_in);
}

double* explicit_simulation() {
    double delta_T = expiration_time/((double)no_of_divisions) ;
    double delta_R = (risk_free_rate - 0.5 * pow(volatility, 2)) * delta_T;
    double delta_SD = volatility * sqrt(delta_T);

    // by sharing random variables we create 4 paths

    
    R = exp(-risk_free_rate * expiration_time);
    double* average=new double[4];
    average[0] = 0;
    average[1] = 0;
    average[2] = 0;
    average[3] = 0;
    for (int i = 0; i < (no_of_trials); i++) {
        int signal1 = 0;
        int signal2 = 0;
        int signal3 = 0;
        int signal4 = 0;
        double current_stock_price1 = initial_stock_price;
        double current_stock_price2 = initial_stock_price;
        double current_stock_price3 = initial_stock_price;
        double current_stock_price4 = initial_stock_price;
        for (int j = 0; j < no_of_divisions; j++) {
            // create the unit normal variates using the Box-Muller Transform
            double x = get_uniform();
            double y = get_uniform();
            double a = sqrt(-2.0 * log(x)) * cos(6.283185307999998 * y);
            double b = sqrt(-2.0 * log(x)) * sin(6.283185307999998 * y);

            current_stock_price1 = current_stock_price1* exp(delta_R + delta_SD * a);
            current_stock_price2 = current_stock_price2 * exp(delta_R - delta_SD * a);
            current_stock_price3 = current_stock_price3 * exp(delta_R + delta_SD * b);
            current_stock_price4 = current_stock_price4 * exp(delta_R - delta_SD * b);
            
            if (current_stock_price1 <= barrier_price) signal1 = 1;
            if (current_stock_price2 <= barrier_price) signal2 = 1;
            if (current_stock_price3 <= barrier_price) signal3 = 1;
            if (current_stock_price4 <= barrier_price) signal4 = 1;
        }
        
        
        double temp1 = current_stock_price1;
        double temp2 = current_stock_price2;
        double temp3 = current_stock_price3;
        double temp4 = current_stock_price4;
        double put1 = current_stock_price1;
        double put2 = current_stock_price2;
        double put3 = current_stock_price3;
        double put4 = current_stock_price4;
        double put11 = current_stock_price1;
        double put21 = current_stock_price2;
        double put31 = current_stock_price3;
        double put41 = current_stock_price4;

        if (signal1 == 1) {
            current_stock_price1 = 0;
            put1 = 0;
        }
        else {
            current_stock_price1 = max((double)0, (current_stock_price1 - strike_price));
            put1= max((double)0, (strike_price-put1));
        }
        
        if (signal2 == 1) {
            current_stock_price2 = 0;
            put2 = 0;
        }
        else {
            current_stock_price2 = max((double)0, (current_stock_price2 - strike_price));
            put2 = max((double)0, (strike_price - put2));
        }
        
        if (signal3 == 1) {
            current_stock_price3 = 0;
            put3 = 0;
        }
        else {
            current_stock_price3 = max((double)0, (current_stock_price3 - strike_price));
            put3 = max((double)0, (strike_price - put3));
        }
        
        if (signal4 == 1) {
            current_stock_price4 = 0;
            put4 = 0;
        }
        else {
            current_stock_price4 = max((double)0, (current_stock_price4 - strike_price));
            put4 = max((double)0, (strike_price - put4));
        }

        
        if (temp1 <= barrier_price) {
            temp1 = 0;
        }
        else {
            temp1 = max((double)0, (temp1 - strike_price)) * (1 - exp((-2.0 / (volatility * volatility * expiration_time)) *
                log(initial_stock_price / barrier_price) *
                log(temp1 / barrier_price)));
        }
        
        if (temp2 <= barrier_price) {
            temp2 = 0;
        }
        else {
            temp2 = max((double)0, (temp2 - strike_price)) * (1 - exp((-2.0 / (volatility * volatility * expiration_time)) *
                log(initial_stock_price / barrier_price) *
                log(temp2 / barrier_price)));
        }
        
        if (temp3 <= barrier_price) {
            temp3 = 0;
        }
        else {
            temp3 = max((double)0, (temp3 - strike_price)) * (1 - exp((-2.0 / (volatility * volatility * expiration_time)) *
                log(initial_stock_price / barrier_price) *
                log(temp3 / barrier_price)));
        }
        
        if (temp4 <= barrier_price) {
            temp4 = 0;
        }
        else {
            temp4 = max((double)0, (temp4 - strike_price)) * (1 - exp((-2.0 / (volatility * volatility * expiration_time)) *
                log(initial_stock_price / barrier_price) *
                log(temp4 / barrier_price)));
        }
        

        if (put11 <= barrier_price) {
            put11 = 0;
        }
        else {
            put11 = max((double)0, (strike_price- put11)) * (1 - exp((-2.0 / (volatility * volatility * expiration_time)) *
                log(initial_stock_price / barrier_price) *
                log(put11 / barrier_price)));
            

        }
        
        if (put21 <= barrier_price) {
            put21 = 0;
        }
        else {
            put21 = max((double)0, (strike_price - put21)) * (1 - exp((-2.0 / (volatility * volatility * expiration_time)) *
                log(initial_stock_price / barrier_price) *
                log(put21 / barrier_price)));
        }
        
        if (put31 <= barrier_price) {
            put31 = 0;
        }
        else {
            put31 = max((double)0, (strike_price - put31)) * (1 - exp((-2.0 / (volatility * volatility * expiration_time)) *
                log(initial_stock_price / barrier_price) *
                log(put31 / barrier_price)));
        }
        
        if (put41 <= barrier_price) {
            put41 = 0;
        }
        else {
            put41 = max((double)0, (strike_price - put41)) * (1 - exp((-2.0 / (volatility * volatility * expiration_time)) *
                log(initial_stock_price / barrier_price) *
                log(put41 / barrier_price)));
        }
        
        
        average[0] = average[0] + (current_stock_price1+ current_stock_price2+ current_stock_price3+ current_stock_price4)/4;
        average[1] = average[1] + (temp1+temp2+temp3+temp4)/4;
        average[2] = average[2] + (put1 + put2 + put3 + put4)/4;
        average[3] = average[3] + (put11 + put21 + put31 + put41)/4;

    }
    average[0] = R * average[0] / ((double)no_of_trials);
    average[1] = R * average[1] / ((double)no_of_trials);
    average[2] = R * average[2] / ((double)no_of_trials);
    average[3] = R * average[3] / ((double)no_of_trials);

    return average;

}


int main(int argc, char* argv[]) {
    sscanf(argv[1], "%lf", &expiration_time);
    sscanf(argv[2], "%lf", &risk_free_rate);
    sscanf(argv[3], "%lf", &volatility);
    sscanf(argv[4], "%lf", &initial_stock_price);
    sscanf(argv[5], "%lf", &strike_price);
    sscanf(argv[6], "%d", &no_of_trials);
    sscanf(argv[7], "%d", &no_of_divisions);
    sscanf(argv[8], "%lf", &barrier_price);
    cout << "----------------------------------" << endl;
    cout << "European Down-and-Out Continuous Barrier Options Pricing via Monte Carlo Simulation" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "Barrier Price = " << barrier_price << endl;
    cout << "Number of Trials = " << no_of_trials << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "----------------------------------" << endl;
    cout << "----------------------------------" << endl;
    
    double* price=explicit_simulation();
    
    cout << "The average Call Price by explicit simulation = " << price[0] << endl;
    cout << "The call price using the (1-p)-adjustment term = " << price[1] << endl;
    cout << "Theoretical Call Price = " << closed_form_down_and_out_european_call_option() << endl;
    cout << "----------------------------------" << endl;
    cout << "----------------------------------" << endl;
    cout << "The average Put Price by explicit simulation = " << price[2] << endl;
    cout << "The put price using the (1-p)-adjustment term = " << price[3] << endl;
    cout << "Theoretical Put Price = " << closed_form_down_and_out_european_put_option() << endl;
    cout << "----------------------------------" << endl;
    return 0;
}

