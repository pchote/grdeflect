/*
 * Copyright 2009, 2012 Paul Chote
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * To link against pgplot, need to complile in two steps:
 * gcc -c main.c
 * gfortran main.o -o test -lm -L/usr/X11R6/lib -lX11 -L/sw/lib -laquaterm -Wl,-framework -Wl,Foundation -L/sw/lib/pgplot -lcpgplot -lpgplot -lpng -lcurses
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <ncurses.h>
#include <cpgplot.h>

const char *versionString = "v1.1 :: 16/03/2012";

#define PI 3.14159265
const double aspect = 1.3;

double u;
double du;
double phi;
double dphi;

/*
 * Functions for the 4th-order Runga-kutta algorithm
 */

// Derivative of u in terms of phi, u, du
double diffu(double phi, double u, double du)
{
    return du;
}

// Derivative of du in terms of phi, u, du
double diffdu(double phi, double u, double du)
{
    // d^2u/dphi^2 + u = 3/2*Rg*u^2
    // if u is dimensionless (u = u*Rg), the Rg drops out
    return (3.0/2)*u*u-u;
}


// Shift position by an amount dphi using 4th-order Runge-Kutta algorithm
void step()
{
    // Evaluate derivatives at the start point
    double k1u = dphi * diffu(phi, u, du);
    double k1du = dphi * diffdu(phi, u, du);

    // Calculate values at intermediate points
    double k2u = dphi * diffu(phi + dphi/2.0, u + k1u/2.0, du + k1du/2.0);
    double k2du = dphi * diffdu(phi + dphi/2.0, u + k1u/2.0, du + k1du/2.0);

    double k3u = dphi * diffu(phi + dphi/2.0, u + k2u/2.0, du + k2du/2.0);
    double k3du = dphi * diffdu(phi + dphi/2.0, u + k2u/2.0, du + k2du/2.0);

    double k4u = dphi * diffu(phi + dphi/2.0, u + k3u/2.0, du + k3du/2.0);
    double k4du = dphi * diffdu(phi + dphi/2.0, u + k3u/2.0, du + k3du/2.0);

    // guess values at end point
    phi += dphi;
    u   += (k1u + 2.0 * k2u + 2.0 * k3u + k4u)/6.0;
    du  += (k1du + 2.0 * k2du + 2.0 * k3du + k4du)/6.0;
}

/*
 * Fire a photon past a mass with an initial separatation of yi on the y axis
 */
int raytrace(double yi, float zoom, bool saveToFile)
{
    // Create an array to hold the calculated data points
    // Start with 100,000 entries which should be sufficient for most cases.
    // Array will expand itself if there is too many entries (will happen if the path comes
    // close to the event horizon and orbits the mass a few times)
    int dataArraySize = 100000;
    double *x = malloc(dataArraySize*sizeof(double));
    double *y = malloc(dataArraySize*sizeof(double));
    double *bigX, *bigY;
    int curEntry = 0;

    // Start the photon at least 1000Rs away;
    double startX = (zoom*aspect > 1000) ? zoom*aspect : 1000;
    phi = PI - atan2(yi, startX);

    // Initial conditions
    u   = sin(phi)/yi;
    du  = cos(phi)/yi;
    dphi = -(0.001/180)*PI; // decrease phi by 0.001 degrees each step

    // Record the closest separation and angle
    double bsq = yi*yi;
    double bphi = PI/2;

    // Calculate path
    bool failed = false;
    bool hires = false;
    while (true)
    {
        // Path collapsing to center - stop calculating before the numbers blow up
        if (zoom*u > 1000)
        {
            failed = true;
            break;
        }
        x[curEntry] = cos(phi)/u;
        y[curEntry] = sin(phi)/u;

        // Squared distance to the center of mass
        double dsq = x[curEntry]*x[curEntry] + y[curEntry]*y[curEntry];
        if (dsq < bsq)
        {
            bsq = dsq;
            bphi = phi;
        }

        // Increase accuracy when d < 2 b
        if (!hires && dsq < 4*yi*yi)
        {
            dphi /= 25;
            hires = true;
        }
        else if (hires && dsq > 4*yi*yi)
        {
            dphi *= 25;
            hires = false;
        }

        // Path has (at least) reached viewport position, and has passed outside the viewport
        // Calculate points to at least 10rs to ensure the path is straight so provide an accurate estimate of deflection angle
        // Also ensure we have at least one point outside the viewport for the line to finish on
        if (phi < PI/2 && (dsq > 100) && ((abs(x[curEntry-1]) > zoom*aspect || abs(y[curEntry-1]) > zoom)))
            break;

        curEntry++;

        // Too many entries, need to expand the data array
        if (curEntry == dataArraySize)
        {
            //printf("Expanding data array %d -> %d entries\n", dataArraySize, 10*dataArraySize);
            bigX = malloc(10*dataArraySize*sizeof(double));
            bigY = malloc(10*dataArraySize*sizeof(double));

            int i;
            for (i = 0; i < dataArraySize; i++)
            {
                bigX[i] = x[i];
                bigY[i] = y[i];
            }
            dataArraySize *= 10;
            free(x);
            free(y);
            x = bigX;
            y = bigY;
        }

        step();
    }

    // Closest separation between photon and center of mass
    double b = sqrt(bsq);
    double deflectAngle = 0;
    if (!failed)
    {
        // Estimate the path as a kinked line, with the kink
        // at the intersection of the line along bphi with yi
        double xi = yi*tan(PI/2 - bphi);

        // For b > 3 estimate deflection from kinked line,
        // otherwise estimate from tangent at closest passing
        deflectAngle = (b > 3) ? -atan2(y[curEntry-1] - yi, x[curEntry-1] - xi) : 2*(PI/2-bphi);
        // Convert to degrees
        deflectAngle *= 180/PI;

        //printf("deflection: %f deg\n", deflectAngle);
        //double bkink = sqrt(xi*xi+yi*yi);
        //printf("theoretical: %f theoretical kink: %f\n", 2/b*180/PI, 2/bkink*180/PI);
        //printf("\ndifference: %4.2f%%\n", (deflectAngle - 2/b*180/PI)/deflectAngle*100);
    }

    // Draw path to screen or ps file
    char filename[20];
    char device[25];

    if (saveToFile)
    {
        sprintf(filename, "path-%.3f.ps",yi);
        sprintf(device, "%s/PS",filename);
    }
    else
        strcpy(device, "9/xs");

    if(cpgbeg(0, device, 1, 1) != 1)
        return EXIT_FAILURE;

    // Setup page, draw axes
    cpgpage();
    cpgsvp(0.05, 0.95, 0.1, 0.95);
    cpgwnad(-aspect*zoom, aspect*zoom, -zoom, zoom);
    cpgbox("bcnts", 0.0, 0, "bcntsv", 0.0, 0);

    cpgscf(2); // Set font style to roman
    cpgslw(2); // Set line width 2

    // Plot title
    char *titleString;
    asprintf(&titleString, "Path for a photon with impact parameter %.3f", yi);
    cpgmtxt("t", 1.0, 0.5, 0.5, titleString);
    free(titleString);

    // Plot event horizon
    cpgsfs(2); // Set fill style to outline
    cpgsls(2); // Set line style to dashed
    cpgsci(2); // Set colour to red
    cpgcirc(0, 0, 1);

    // Create a float array for pgplot
    float *xf = malloc(curEntry*sizeof(float));
    float *yf = malloc(curEntry*sizeof(float));
    for (int i = 0; i < curEntry; i++)
    {
        xf[i] = x[i];
        yf[i] = y[i];
    }

    cpgsci(1); // White
    char *deflabel, *blabel;
    if (failed)
        asprintf(&deflabel, "Total deflection: Captured!");
    else if (deflectAngle > 0.1)
        asprintf(&deflabel, "Total deflection: \\gD\\gh = %3.3f deg", deflectAngle);
    else
        asprintf(&deflabel, "Total deflection: \\gD\\gh = %1.3e deg", deflectAngle);

    if (failed)
        asprintf(&blabel, "Closest approach: Captured!");
    else
        asprintf(&blabel, "Closest approach: b = %0.3f", b);

    cpgmtxt("b", 3, 0, 0, blabel);
    cpgmtxt("b", 3, 1, 1, deflabel);

    free(blabel);
    free(deflabel);

    // Plot path
    cpgsls(1); // Set line style to normal
    cpgsci(7); // Set colour to yellow
    cpgline(curEntry,xf,yf);
    cpgend();

    if (saveToFile)
    {
        clear();
        printw("Saved to file '%s'\n",filename);
        printw("\n(press any key to continue)\n");
        getch();
    }

    free(x);
    free(y);
    return EXIT_SUCCESS;
}

/*
 * Prints a nicely formatted command list
 */
void printHelp()
{
    clear();
    attron(A_BOLD);
    printw("grdeflect ");
    attroff(A_BOLD);
    printw(versionString);
    printw("\nby Paul Chote for Victoria University of Wellington\n");
    int row = 3;
    int cola = 3;
    int colb = 30;

    move(row,0);
    attron(A_BOLD);
    printw("OVERVIEW:");
    attroff(A_BOLD);
    move(++row,cola);
    printw("Simulates the path of a photon travelling past a gravitational body.");
    move(++row,cola);
    printw("A photon is fired horizontally from the left of the screen, with an");
    move(++row,cola);
    printw("impact parameter specified by the FIRE command (given as a ");
    move(++row,cola);
    printw("multiple of the Schwarzschild radius, Rs).");

    row += 2;
    move(row,cola);
    printw("The path of the photon will be bent, or it will be captured if it");
    move(++row,cola);
    printw("passes within 1Rs of the center of mass.");

    row += 2;
    move(row,cola);
    printw("Most objects have a physical size much larger than their Schwarzschild");
    move(++row,cola);
    printw("radius, which prevents photons from passing this close. The exception");
    move(++row,cola);
    printw("are Black Holes, where 1Rs defines their event horizon.");

    row += 2;

    move(row,0);
    attron(A_BOLD);
    printw("COMMANDS:");
    attroff(A_BOLD);

    move(++row,cola);
    attron(A_BOLD);
    printw("?");
    attroff(A_BOLD);
    move(row,colb);
    printw("Display this description.");

    move(++row,cola);
    attron(A_BOLD);
    printw("FIRE <impact parameter>");
    attroff(A_BOLD);
    move(row,colb);
    printw("Fire a photo with specified impact parameter.");

    move(++row,cola);
    attron(A_BOLD);
    printw("VIEWPORT <maximum radius>");
    attroff(A_BOLD);
    move(row,colb);
    printw("Adjust the zoom of the viewport.");

    move(++row,cola);
    attron(A_BOLD);
    printw("PRINT");
    attroff(A_BOLD);
    move(row,colb);
    printw("Save current graph to file.");

    move(++row,cola);
    attron(A_BOLD);
    printw("QUIT");
    attroff(A_BOLD);
    move(row,colb);
    printw("Exit the program.");

    move(row+2, 0);
    printw("(press any key to continue)");

    refresh();

    getch();

    clear();
}

void printError(char *msg)
{
    clear();
    printw(msg);
    printw("\n\n(press any key to continue)");

    getch();
}

/*
 * I am using the ncurses library to improve the user experience, at the cost
 * of slightly increasing the complexity of the code.
 *
 * ncurses lets the program "take over" the terminal window so we can control what is shown,
 * instead of having a history of commands shown, and lets us use formatting etc in the text.
 * The main difference is that we must use printw instead of printf to write text to the screen
 *
 */

int main(int argc, char *argv[])
{
    initscr();              /* start the curses mode */
    echo();                 /* show user input */

    double b = 3.0f;
    float zoom = 5.0f;

    double inputb;
    float inputzoom;
    bool doRaytrace = false;
    char inputBuffer[20];

    while (true)
    {
        if (doRaytrace)
        {
            raytrace(b, zoom, false);
            doRaytrace = 0;
        }
        clear();            /* clear the screen */
        printw("Enter Command (type \"?\" for help):\n");

        getstr(inputBuffer); // Get the user's reply

        if (strcasecmp("QUIT", inputBuffer) == 0)
        {
            break;
        }

        if (strcasecmp("?", inputBuffer) == 0)
        {
            printHelp();
            continue;
        }

        if (strcasecmp("PRINT", inputBuffer) == 0)
        {
            raytrace(b, zoom, true);
            continue;
        }

        if (strncasecmp("FIRE ", inputBuffer,5) == 0)
        {
            inputb = -1;
            sscanf(inputBuffer,"%*s %lf",&inputb);
            if (b != inputb)
            {
                if (inputb >= 0)
                {
                    b = inputb;
                    doRaytrace = 1;
                }
                else
                    printError("Invalid impact parameter specified.");
            }
            continue;
        }

        if (strncasecmp("VIEWPORT ", inputBuffer,9) == 0)
        {
            inputzoom = -1;
            sscanf(inputBuffer,"%*s %f",&inputzoom);
            if (zoom != inputzoom)
            {
                if (inputzoom > 0)
                {
                    zoom = inputzoom;
                    doRaytrace = 1;
                }
                else
                    printError("Invalid zoom specified.");
            }
            continue;
        }
        printError("Unknown command");
    }

    endwin();               /* end curses mode */
    return EXIT_SUCCESS;
}
