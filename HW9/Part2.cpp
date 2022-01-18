//
//  Part2.cpp
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

double* explicit_simulation() {
    double delta_T = expiration_time / ((double)no_of_sampling_instants);
    double delta_R = (risk_free_rate - 0.5 * pow(volatility, 2)) * delta_T;
    double delta_SD = volatility * sqrt(delta_T);

    // by sharing random variables we create 4 paths


    R = exp(-risk_free_rate * expiration_time);
    double* average = new double[4];
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
        for (int j = 0; j < no_of_sampling_instants; j++) {
            // create the unit normal variates using the Box-Muller Transform
            double x = get_uniform();
            double y = get_uniform();
            double a = sqrt(-2.0 * log(x)) * cos(6.283185307999998 * y);
            double b = sqrt(-2.0 * log(x)) * sin(6.283185307999998 * y);
    
            current_stock_price1 = current_stock_price1 * exp(delta_R + delta_SD * a);
            current_stock_price2 = current_stock_price2 * exp(delta_R - delta_SD * a);
            current_stock_price3 = current_stock_price3 * exp(delta_R + delta_SD * b);
            current_stock_price4 = current_stock_price4 * exp(delta_R - delta_SD * b);

            if (current_stock_price1 <= barrier_price)
                signal1 = 1;
            if (current_stock_price2 <= barrier_price)
                signal2 = 1;
            if (current_stock_price3 <= barrier_price)
                signal3 = 1;
            if (current_stock_price4 <= barrier_price)
                signal4 = 1;
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
            put1 = max((double)0, (strike_price - put1));
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
        
        
        double* mean=new double[no_of_sampling_instants-1];
        double* variance = new double[no_of_sampling_instants-1];
        double p_not_hit;

        
        if (temp1 <= barrier_price) {
            temp1 = 0;
        }
        else {
            p_not_hit = 1.0;
            for (int j = 0; j < no_of_sampling_instants-1; j++) {
                mean[j] = initial_stock_price +
                    ((float)j+1) / ((float)no_of_sampling_instants) * (temp1 - initial_stock_price);
                variance[j] = (((float)j+1) / ((float)no_of_sampling_instants)) * expiration_time *
                    (1.0 - ((float)(j+1)) / ((float)no_of_sampling_instants));
            }
            for (int j = 0; j < no_of_sampling_instants-1; j++)
                p_not_hit *=
                (1.0 - N((barrier_price - mean[j]) / sqrt(variance[j])));

            temp1 = max((double)0, (temp1 - strike_price)) * p_not_hit;
        }
        
        if (temp2 <= barrier_price) {
            temp2 = 0;
        }
        else {
            p_not_hit = 1.0;
            for (int j = 0; j < no_of_sampling_instants-1; j++) {
                mean[j] = initial_stock_price +
                    ((float)j+1) / ((float)no_of_sampling_instants) * (temp2 - initial_stock_price);
                variance[j] = (((float)j+1) / ((float)no_of_sampling_instants)) * expiration_time *
                    (1.0 - ((float)(j+1)) / ((float)no_of_sampling_instants));
            }
            for (int j = 0; j < no_of_sampling_instants-1; j++)
                p_not_hit *=
                (1.0 - N((barrier_price - mean[j]) / sqrt(variance[j])));

            temp2 = max((double)0, (temp2 - strike_price)) * p_not_hit;
        }
        
        if (temp3 <= barrier_price) {
            temp3 = 0;
        }
        else {
            p_not_hit = 1.0;
            for (int j = 0; j < no_of_sampling_instants-1; j++) {
                mean[j] = initial_stock_price +
                    ((float)j+1) / ((float)no_of_sampling_instants) * (temp3 - initial_stock_price);
                variance[j] = (((float)j+1) / ((float)no_of_sampling_instants)) * expiration_time *
                    (1.0 - ((float)(j+1)) / ((float)no_of_sampling_instants));
            }
            for (int j = 0; j < no_of_sampling_instants-1; j++)
                p_not_hit *=
                (1.0 - N((barrier_price - mean[j]) / sqrt(variance[j])));

            temp3 = max((double)0, (temp3 - strike_price)) * p_not_hit;
        }
        
        if (temp4 <= barrier_price) {
            temp4 = 0;
        }
        else {
            p_not_hit = 1.0;
            for (int j = 0; j < no_of_sampling_instants-1; j++) {
                mean[j] = initial_stock_price +
                    ((float)j+1) / ((float)no_of_sampling_instants) * (temp4 - initial_stock_price);
                variance[j] = (((float)j+1) / ((float)no_of_sampling_instants)) * expiration_time *
                    (1.0 - ((float)(j+1)) / ((float)no_of_sampling_instants));
            }
            for (int j = 0; j < no_of_sampling_instants-1; j++)
                p_not_hit *=
                (1.0 - N((barrier_price - mean[j]) / sqrt(variance[j])));

            temp4 = max((double)0, (temp4 - strike_price)) * p_not_hit;
        }

        if (put11 <= barrier_price) {
            put11 = 0;
        }
        else {
            p_not_hit = 1.0;
            for (int j = 0; j < no_of_sampling_instants-1; j++) {
                mean[j] = initial_stock_price +
                    ((float)j+1) / ((float)no_of_sampling_instants) * (put11 - initial_stock_price);
                variance[j] = (((float)j+1) / ((float)no_of_sampling_instants)) * expiration_time *
                    (1.0 - ((float)(j+1)) / ((float)no_of_sampling_instants));
            }
            for (int j = 0; j < no_of_sampling_instants-1; j++)
                p_not_hit *=
                (1.0 - N((barrier_price - mean[j]) / sqrt(variance[j])));

            put11 = max((double)0, (strike_price-put11)) * p_not_hit;
        }
        
        if (put21 <= barrier_price) {
            put21 = 0;
        }
        else {
            p_not_hit = 1.0;
            for (int j = 0; j < no_of_sampling_instants-1; j++) {
                mean[j] = initial_stock_price +
                    ((float)j+1) / ((float)no_of_sampling_instants) * (put21 - initial_stock_price);
                variance[j] = (((float)j+1) / ((float)no_of_sampling_instants)) * expiration_time *
                    (1.0 - ((float)(j+1)) / ((float)no_of_sampling_instants));
            }
            for (int j = 0; j < no_of_sampling_instants-1; j++)
                p_not_hit *=
                (1.0 - N((barrier_price - mean[j]) / sqrt(variance[j])));

            put21 = max((double)0, (strike_price - put21)) * p_not_hit;
        }
        
        if (put31 <= barrier_price) {
            put31 = 0;
        }
        else {
            p_not_hit = 1.0;
            for (int j = 0; j < no_of_sampling_instants-1; j++) {
                mean[j] = initial_stock_price +
                    ((float)j+1) / ((float)no_of_sampling_instants) * (put31 - initial_stock_price);
                variance[j] = (((float)j+1) / ((float)no_of_sampling_instants)) * expiration_time *
                    (1.0 - ((float)(j+1)) / ((float)no_of_sampling_instants));
            }
            for (int j = 0; j < no_of_sampling_instants-1; j++)
                p_not_hit *=
                (1.0 - N((barrier_price - mean[j]) / sqrt(variance[j])));

            put31 = max((double)0, (strike_price - put31)) * p_not_hit;
        }
        
        if (put41 <= barrier_price) {
            put41 = 0;
        }
        else {
            p_not_hit = 1.0;
            for (int j = 0; j < no_of_sampling_instants-1; j++) {
                mean[j] = initial_stock_price +
                    ((float)j+1) / ((float)no_of_sampling_instants) * (put41 - initial_stock_price);
                variance[j] = (((float)j+1) / ((float)no_of_sampling_instants)) * expiration_time *
                    (1.0 - ((float)(j+1)) / ((float)no_of_sampling_instants));
            }
            for (int j = 0; j < no_of_sampling_instants-1; j++)
                p_not_hit *=
                (1.0 - N((barrier_price - mean[j]) / sqrt(variance[j])));

            put41 = max((double)0, (strike_price - put41)) * p_not_hit;
        }
        average[0] = average[0] + (current_stock_price1 + current_stock_price2 + current_stock_price3 + current_stock_price4)/4;
        average[1] = average[1] + (temp1 + temp2 + temp3 + temp4)/4;
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
    sscanf(argv[7], "%d", &no_of_sampling_instants);
    sscanf(argv[8], "%lf", &barrier_price);
    cout << "----------------------------------" << endl;
    cout << "European Down-and-Out Discrete Barrier Options Pricing via Monte Carlo Simulation" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "Barrier Price = " << barrier_price << endl;
    cout << "Number of Trials = " << no_of_trials << endl;
    cout << "Number of Discrete Barriers = " << no_of_sampling_instants << endl;
    cout << "----------------------------------" << endl;
    
    double* price = explicit_simulation();

    cout << "The average Call Price by explicit simulation of the price paths = " << price[0] << endl;
    cout << "The average Call Price with Brownian-Bridge correction on the final price = " << price[1] << endl;
    cout << "The average Put Price by explicit simulation of the price paths = " << price[2] << endl;
    cout << "The average Put Price with Brownian-Bridge correction on the final price = " << price[3] << endl;
    cout << "----------------------------------" << endl;
    return 0;
}
