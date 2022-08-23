//
// Created by Konstantinos Dalampekis on 17/8/22.
//

/*
 * sx0 -> the initial value of fibre length.
 * In this code k indicates the index of the differential equations.
 * k = 1 -> Signal, k = 2 -> Pa, k = 3 -> Pn.
 * During the input you will be asked to enter the number of steps.
 * The steps are the points where the output will show the values along the length.
 * The output is standard to show only the values for all wavelengths at x = xfinal.
 * So, if you want to see also the values along the length you can choose a step size.
 * Example if steps = 1 -> it will show the values of Power for each wavelength and x = 0,1,2,3...xfinal
 * if steps = 0.1 -> x = 0.1,0.2,...xfinal.
 * if you don't care about the values along the fibre length set steps = 0.
 * Then, if you choose k = 3, the input will ask you to set values for parameters c, b.
 * c and b are β and γ at the double cladding equation.
 * the output will show you the values of each Power as well as,
 * the maximum value of Signal
 * the minimum value of Pa
 * the minimum value of Pn
 * the minimum value of Ptotal
 * across all the wavelengths and fibre length.
 */

#include <iostream>
#include <cmath>
#include "fstream"
#include "string"

using namespace std;

double  Gk[3]; // overlap integral
int k;
double sa[1000], se[1000]; /* absorption and emission sections */
int l, lin = 1, lfin = 362;
double wav[1000];// wavelength //
double H = 6.626e-34, z = 4.2*pow(10,15), uk=1, lk = 0.00; // Ni = E[ni]
double nt = 1.00e+25; // Yterbioum ions concentration
long double y[1000][1000]; // y = y[k][z]


double vk(int wav0) /* frequency */
{
    return (3e+8)/(wav[wav0]*(1e-9));
}
double ak(int k1, int wav1) /* loss specrtum */
{
    return sa[wav1]*Gk[k1]*nt;
}
double gk(int k2, int wav2) /* gain spectrum */
{
    return se[wav2]*Gk[k2]*nt;
}

int main() {

    int maxSl = 1, maxP2l = 1, maxP3l = 1, maxP4l = 1; /* wavelength -> max Signal, max non-overlapPump and max overlap Pump */
    double c, bb;
    int x0, xn, steps2;
    double steps;

    /* copying the values of wavelength from text file */
    short loop = 1;
    string line;
    ifstream myfile("/Volumes/Mac Os - Δεδομένα/Clion_Projects/Fibre_Laser/Wavelengths.txt");
    if (myfile.is_open()) {
        while (!myfile.eof()) {
            getline(myfile, line); /* reading the values as a string  */
            double num = stod(line); /* convert the string to double */
            wav[loop] = num;
            loop++;
        }
        myfile.close();
    }

    /* copying the values of absorption section from text file */
    loop = 1;
    ifstream myfile1("/Volumes/Mac Os - Δεδομένα/Clion_Projects/Fibre_Laser/Absorption_section.txt");
    if (myfile1.is_open()) {
        while (!myfile1.eof()) {
            getline(myfile1, line); /* reading the values as a string  */
            double num = stod(line); /* convert the string to double */
            sa[loop] = num;
            loop++;
        }
        myfile1.close();
    }

    /* copying the values of emission section from text file */
    loop = 1;
    ifstream myfile2("/Volumes/Mac Os - Δεδομένα/Clion_Projects/Fibre_Laser/Emission_section.txt");
    if (myfile2.is_open()) {
        while (!myfile2.eof()) {
            getline(myfile2, line); /* reading the values as a string  */
            double num = stod(line); /* convert the string to double */
            se[loop] = num;
            loop++;
        }
        myfile2.close();
    }

    /* Entering the initial conditions and constants */
    cout << "Enter Initial Condition" << endl;
    cout << "x0 = ";
    cin >> x0;
    int initial = x0;
    cout<<"If you want to use the amplifier mode set k = 2"<< endl;
    cout<<"if you want to use the double cladding mode set k = 3"<< endl;
    cout<<"Enter value of coefficient k: ";
    cin >> k;
    if (k <= 1) {
        cout << "Enter a minimum value 2, try again: ";
        cin >> k;}
    y[0][0] = 0;
    cout << "Enter the value of Signal when z = "; cout << x0 << ":" << endl;
    cin >> y[1][x0];
    cout << "Enter the value of Pa when z = "; cout << x0;
    cout << " and wavelength = " << wav[lin] << ":" << endl;
    cin >> y[2][x0];
    if (k == 3){
        cout << "Enter the value of Pn when z = "; cout << x0;
        cout << " and wavelength = " << wav[lin] << ":" << endl;
        cin >> y[3][x0];
    }
    cout << "Enter calculation point xn = ";
    cin >> xn;
    if(k == 3) xn = 10*xn;//xn = z
    cout << "if k = 2 enter an integer number of steps" << endl;
    cout << "if k = 3 you can use also real number for steps" << endl;
    cout<<"Enter number of steps(if you do not care about the values in the fiber length then set steps = 0): ";
    if(k == 2){cin >> steps2; }
    if(k == 3){cin >> steps;}
    if(k==3){
        cout<<"Enter value of coefficient c: ";
        cin >> c;
        cout<<"Enter value of coefficient b: ";
        cin >> bb;}

    /* Setting the values of Gk for each k */
    Gk[1] = 0.9;
    Gk[2] = 0.1;

    /* max Signal, max non-overlapPump and max overlap Pump */
    long double max1 = y[k-2][0], max2 = y[k-1][0], max3 = y[k][0];

    /* Runge Kutta Method, calculation and display of result */
    cout << "\nx0\tyn\n";
    cout << "------------------\n";
    if(k==3) {
        /* Calculating step size (h) */
        int hh = 1;
        double h = 0.1;
        int xmax1 = 0, xmax2 = 0, xmax3 = 0, xmax4 = 0;
        int l1 = 280;
        long double aabb = y[k-1][0] + y[k][0];
        long double n;
        for (l = lin; l <= lfin; l++) { /* for loop the wavelength */
            for (int xi = initial; xi <= xn; xi++) { /* for loop the fibre length */
                n = ((y[1][xi] * ak(1, l1) / (H * vk(l1) * z)) + ((y[2][xi] * ak(2, l)) / (H * vk(l) * z))) / (1 +
                                                                                                               ((y[1][xi] *
                                                                                                                 (ak(1,
                                                                                                                     l1) +
                                                                                                                  gk(1,
                                                                                                                     l1)) /
                                                                                                                 (H *
                                                                                                                  vk(l1) *
                                                                                                                  z)) +
                                                                                                                (y[2][xi] *
                                                                                                                 (ak(2, l) +
                                                                                                                  gk(2,
                                                                                                                     l)) /
                                                                                                                 (H * vk(l) *
                                                                                                                  z))));
                y[k - 2][xi + hh] = y[k - 2][xi] + h * (y[k - 2][xi] * (uk * n * (ak(k - 2, l1) + gk(k - 2, l1)) -
                                                                        (uk * (ak(k - 2, l1) + lk))));
                y[k - 1][xi + hh] = y[k - 1][xi] + h * (y[k - 1][xi] *
                                                        (-n * (ak(k - 1, l) + gk(k - 1, l)) - (ak(k - 1, l) + lk) - c) +
                                                        bb * y[k][xi]);
                y[k][xi + hh] = y[k][xi] + h * (c * y[k - 1][xi] - bb * y[k][xi]);

                /* Calculation of max values for Signal and Pumps in function to wavelength */
                if (y[k - 2][xi] > max1 && k >= 1 && xi != 0) {
                    xmax1 = xi;
                    max1 = y[k - 2][xi];
                    maxSl = l;
                }
                if (y[k - 1][xi] < max2 && k >= 2 && xi != 0) {
                    xmax2 = xi;
                    max2 = y[k - 1][xi];
                    maxP2l = l;
                }
                if (y[k][xi] < max3 && k >= 3 && xi != 0) {
                    xmax3 = xi;
                    max3 = y[k][xi];
                    maxP3l = l;
                }
                if (y[k][xi] + y[k - 1][xi] < aabb && xi != 0) {
                    xmax4 = xi;
                    maxP4l = l;
                    aabb = y[k][xi] + y[k - 1][xi];
                }

            }
            /* Displaying the results among the fibre length */
            if(steps!=0) {
                for (int xi = 0; xi <= xn; xi = xi + 10*steps) {
                    cout << xi*0.1  << "\t" << "Signal " << y[k - 2][xi] << " Pa " << y[k - 1][xi] << " Pn "
                         << y[k][xi]
                         << " P "
                         << y[k - 1][xi] + y[k][xi] << endl;
                }
            }
            /* Displaying the results along the wavelength */
            cout << "\nValue of y at x = " << xn * 0.1 << " and at l = " << wav[l] << " is " << "Signal "
                 << y[k - 2][xn]
                 << " Pa " << y[k - 1][xn] << " Pn " << y[k][xn] << " P " << y[k - 1][xn] + y[k][xn] << endl;


        }
        /* Display of max values for each Power */
        if (k >= 1)
            cout << " The Max value of Signal is at l = " << wav[maxSl] << " , x = " << xmax1 * 0.1 << " and Signal = "
                 << max1 << endl;
        if (k >= 2)
            cout << " The Min value of non-overlapping Pump is at l = " << wav[maxP2l] << " , x = " << xmax2 * 0.1
                 << " and Pa = " << max2 << endl;
        if (k >= 3)
            cout << " The Min value of overlapping Pump is at l = " << wav[maxP3l] << " , x = " << xmax3 * 0.1
                 << " and Pn = " << max3 << endl;
        if (k >= 3)
            cout << " The Min value of Sum(Pa,Pn) is at l = " << wav[maxP4l] << " , x = " << xmax4 * 0.1
                 << " and Pa + Pn = " << aabb << endl;
    }
    if( k == 2 ){
        int h = 1;
        long double n1,n2;
        for (int i = initial; i <= xn; i++) { /* for loop the fibre length */
            if (k == 2) { /* Pump and Signal calculation */
                for (int j = 1; j < k; j++) {
                    n1 = ((y[1][i] * ak(1, 280) / (H * vk(280) * z)) + ((y[2][i] * ak(2, 155)) / (H * vk(155) * z))) / (1 +
                                                                                                                        ((y[1][i] *
                                                                                                                          (ak(1,
                                                                                                                              280) +
                                                                                                                           gk(1,
                                                                                                                              280)) /
                                                                                                                          (H *
                                                                                                                           vk(280) *
                                                                                                                           z)) +
                                                                                                                         (y[2][i] *
                                                                                                                          (ak(2,
                                                                                                                              155) +
                                                                                                                           gk(2,
                                                                                                                              155)) /
                                                                                                                          (H *
                                                                                                                           vk(155) *
                                                                                                                           z))));
                    y[k - j][i + h] = y[k - j][i] + (h * ((uk * (ak(1, 280) + gk(1, 280)) * (n1 * y[k - j][i])) -
                                                          (uk * (ak(1, 280) + lk) * y[k - j][i])));
                }
                n2 = ((y[1][i] * ak(1, 280) / (H * vk(280) * z)) + ((y[2][i] * ak(2, 155)) / (H * vk(155) * z))) / (1 +
                                                                                                                    ((y[1][i] *
                                                                                                                      (ak(1,
                                                                                                                          280) +
                                                                                                                       gk(1,
                                                                                                                          280)) /
                                                                                                                      (H *
                                                                                                                       vk(280) *
                                                                                                                       z)) +
                                                                                                                     (y[2][i] *
                                                                                                                      (ak(2,
                                                                                                                          155) +
                                                                                                                       gk(2,
                                                                                                                          155)) /
                                                                                                                      (H *
                                                                                                                       vk(155) *
                                                                                                                       z))));
                y[k][i + h] = y[k][i] + (h * ((uk * (ak(2, 155) + gk(2, 155)) * (n2 * y[k][i])) -
                                              (uk * (ak(2, 155) + lk) * y[k][i])));
            }
            /* Displaying the results along the fibre length */
            if (steps2 != 0) {
                if (i % steps2 == 0 && k == 2)
                    cout << i << "\t" << "Signal " << y[k - 1][i] << "\t" << "Pump " << y[k][i] << endl;
                if (i % steps2 == 0 && k == 1 && steps2 != 0) cout << i << "\t" << "Signal " << y[k][i] << endl;
            }
        }
    }
    return 0;
}