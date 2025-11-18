/*
 * Project Title: The C Optics Calculator: Rings, Gratings, and Rotations
 * Description: A comprehensive tool for calculating optical phenomena including
 * Newton's Rings, Diffraction Gratings, Dispersive Power,
 * Malus's Law, and Specific Rotation.
 *
 * Compilation: gcc optics_calculator.c -o optics -lm
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ================= FUNCTION PROTOTYPES =================

void showMainMenu();
void handleNewtonsRings();
void calculateNewtonsRadius();
void calculateNewtonsWavelength();
void handleDispersivePower();
void handleDiffractionGrating();
void handleMalusLaw();
void handleSpecificRotation();
void printCredits();
double deg2rad(double deg);
double calculateRefractiveIndex(double A_deg, double delta_m_deg);

// ================= MAIN FUNCTION =================

int main() {
    int choice;

    do {
        showMainMenu();
        printf("Enter your choice (0-5): ");
        if (scanf("%d", &choice) != 1) {
            // Handle non-integer input to prevent infinite loops
            while(getchar() != '\n');
            choice = -1;
        }

        switch(choice) {
            case 1:
                handleNewtonsRings();
                break;
            case 2:
                handleDiffractionGrating();
                break;
            case 3:
                handleDispersivePower();
                break;
            case 4:
                handleMalusLaw();
                break;
            case 5:
                handleSpecificRotation();
                break;
            case 0:
                printCredits();
                printf("\nExiting program. Goodbye!\n");
                break;
            default:
                printf("\n[!] Invalid choice. Please try again.\n");
        }

        if (choice != 0) {
            printf("\nPress Enter to return to main menu...");
            while(getchar() != '\n'); // Clear buffer
            getchar(); // Wait for Enter
        }

    } while (choice != 0);

    return 0;
}

// ================= MENU INTERFACE =================

void showMainMenu() {
    // Clear screen command (system dependent, optional)
    // system("cls"); // Windows
    // system("clear"); // Linux/Mac

    printf("\n=======================================================\n");
    printf("   The C Optics Calculator: Rings, Gratings, and Rotations\n");
    printf("=======================================================\n");
    printf(" 1. Newton's Rings (Radius & Wavelength)\n");
    printf(" 2. Diffraction Grating Analysis\n");
    printf(" 3. Dispersive Power of a Prism\n");
    printf(" 4. Malus's Law Simulation\n");
    printf(" 5. Specific Rotation Calculator\n");
    printf(" 0. Exit & Credits\n");
    printf("=======================================================\n");
}

// ================= MODULE 1: NEWTON'S RINGS =================

void handleNewtonsRings() {
    int subChoice;
    printf("\n--- Newton's Rings Module ---\n");
    printf("1. Calculate Radius of nth Ring\n");
    printf("2. Calculate Wavelength from Ring Radius\n");
    printf("Enter choice: ");
    scanf("%d", &subChoice);

    if (subChoice == 1) {
        calculateNewtonsRadius();
    } else if (subChoice == 2) {
        calculateNewtonsWavelength();
    } else {
        printf("Invalid sub-choice.\n");
    }
}

void calculateNewtonsRadius() {
    double lambda, R, n, rn;

    printf("\n[Newton's Ring Radius Calculator]\n");
    printf("Formula: r_n = sqrt(n * lambda * R)\n");

    printf("Enter wavelength of light (in meters, e.g., 589e-9): ");
    scanf("%lf", &lambda);
    printf("Enter ring number (n): ");
    scanf("%lf", &n);
    printf("Enter radius of curvature of lens (R in meters): ");
    scanf("%lf", &R);

    if (n < 0 || R < 0 || lambda < 0) {
        printf("Error: Physical parameters cannot be negative.\n");
        return;
    }

    rn = sqrt(n * lambda * R);
    printf(">> The radius of the %.0f-th ring is: %e meters\n", n, rn);
}

void calculateNewtonsWavelength() {
    double rn, n, R, lambda;

    printf("\n[Newton's Ring Wavelength Calculator]\n");
    printf("Formula: lambda = r_n^2 / (n * R)\n");

    printf("Enter radius of the nth ring (r_n in meters): ");
    scanf("%lf", &rn);
    printf("Enter ring number (n): ");
    scanf("%lf", &n);
    printf("Enter radius of curvature of lens (R in meters): ");
    scanf("%lf", &R);

    if (n <= 0 || R <= 0) {
        printf("Error: n and R must be positive and non-zero.\n");
        return;
    }

    lambda = (rn * rn) / (n * R);
    printf(">> The calculated wavelength of light is: %e meters\n", lambda);
}

// ================= MODULE 2: DIFFRACTION GRATING =================

void handleDiffractionGrating() {
    double theta_deg, theta_rad, lambda_nm, lambda_m, n, d, N, N_per_inch;

    printf("\n--- Diffraction Grating Calculator ---\n");
    printf("Formula: d * sin(theta) = n * lambda\n");

    printf("Enter angle of diffraction (theta in degrees): ");
    scanf("%lf", &theta_deg);
    printf("Enter order of diffraction (n): ");
    scanf("%lf", &n);
    printf("Enter wavelength of light (in nanometers): ");
    scanf("%lf", &lambda_nm);

    // Conversions
    lambda_m = lambda_nm * 1e-9;
    theta_rad = deg2rad(theta_deg);

    // Validation
    if (sin(theta_rad) == 0) {
        printf("Error: sin(theta) is zero. Cannot divide by zero.\n");
        return;
    }
    if (n <= 0) {
        printf("Error: Order (n) must be positive.\n");
        return;
    }

    // Calculations
    d = (n * lambda_m) / sin(theta_rad);
    N = 1.0 / d;
    N_per_inch = N * 0.0254;

    printf("\n=== Results ===\n");
    printf("Calculated grating spacing (d): %e meters\n", d);
    printf("Lines per meter (N):            %e lines/m\n", N);
    printf("Lines per inch:                 %e lines/inch\n", N_per_inch);
}

// ================= MODULE 3: DISPERSIVE POWER =================

double calculateRefractiveIndex(double A_deg, double delta_m_deg) {
    double A = deg2rad(A_deg);
    double delta = deg2rad(delta_m_deg);
    double num = sin((A + delta) / 2.0);
    double den = sin(A / 2.0);

    if (den == 0) return 0.0; // Prevent division by zero
    return num / den;
}

void handleDispersivePower() {
    double twoA, two_delta_r, two_delta_g, two_delta_v;
    double A, delta_r, delta_g, delta_v;
    double mu_r, mu_g, mu_v, omega;

    printf("\n--- Dispersive Power Calculator ---\n");
    printf("Formula: omega = (mu_v - mu_r) / (mu_y - 1)\n");
    printf("Note: Ensure inputs correspond to Red, Green/Yellow, and Violet lines.\n");

    printf("Enter measured 2A (in degrees): ");
    scanf("%lf", &twoA);
    A = twoA / 2.0;

    printf("Enter measured 2*delta_m for RED (in degrees): ");
    scanf("%lf", &two_delta_r);
    printf("Enter measured 2*delta_m for GREEN/YELLOW (in degrees): ");
    scanf("%lf", &two_delta_g);
    printf("Enter measured 2*delta_m for VIOLET (in degrees): ");
    scanf("%lf", &two_delta_v);

    // Divide by 2 to get actual angle of minimum deviation
    delta_r = two_delta_r / 2.0;
    delta_g = two_delta_g / 2.0;
    delta_v = two_delta_v / 2.0;

    mu_r = calculateRefractiveIndex(A, delta_r);
    mu_g = calculateRefractiveIndex(A, delta_g);
    mu_v = calculateRefractiveIndex(A, delta_v);

    if (mu_g <= 1.0) {
        printf("Error: Refractive index of mean ray must be > 1.\n");
        return;
    }

    omega = (mu_v - mu_r) / (mu_g - 1.0);

    printf("\n=== Results ===\n");
    printf("Prism Angle (A):          %.2f degrees\n", A);
    printf("Refractive Index (Red):    %.6f\n", mu_r);
    printf("Refractive Index (Mean):   %.6f\n", mu_g);
    printf("Refractive Index (Violet): %.6f\n", mu_v);
    printf("Dispersive Power (omega):  %e\n", omega);
}

// ================= MODULE 4: MALUS'S LAW =================

void handleMalusLaw() {
    double I0, angle, intensity;
    int step;

    printf("\n--- Malus's Law Simulation ---\n");
    printf("Formula: I = I0 * cos^2(theta)\n");

    printf("Enter initial intensity (I0): ");
    scanf("%lf", &I0);
    printf("Enter step size for table (in degrees, e.g., 30): ");
    scanf("%d", &step);

    if (step <= 0) {
        printf("Error: Step size must be positive.\n");
        return;
    }

    printf("\n Angle (deg) | Intensity\n");
    printf("-------------|-----------\n");

    for (angle = 0; angle <= 180; angle += step) {
        double rad = deg2rad(angle);
        intensity = I0 * pow(cos(rad), 2);
        printf(" %6.1f      | %.4f\n", angle, intensity);
    }
}

// ================= MODULE 5: SPECIFIC ROTATION =================

void handleSpecificRotation() {
    double theta, L, c, s;

    printf("\n--- Specific Rotation Calculator ---\n");
    printf("Formula: S = theta / (L * c)\n");

    printf("Enter observed rotation (theta in degrees): ");
    scanf("%lf", &theta);
    printf("Enter length of tube (L in decimeters/dm): ");
    scanf("%lf", &L);
    printf("Enter concentration (c in g/mL): ");
    scanf("%lf", &c);

    if (L <= 0 || c <= 0) {
        printf("Error: Tube length and concentration must be positive values.\n");
        return;
    }

    s = theta / (L * c);
    printf(">> Specific Rotation [alpha]: %.3f degrees dm^-1 (g/mL)^-1\n", s);
}

// ================= UTILITIES & CREDITS =================

double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}

void printCredits() {
    printf("\n=======================================================\n");
    printf("                    CREDITS & REFERENCES               \n");
    printf("=======================================================\n");
    printf("Student Contributors:\n");
    printf("1. NEWTON'S RING:           Rudra (992501030014)\n");
    printf("2. DISPERSIVE POWER:        Aanya (992501030005)\n");
    printf("3. DIFFRACTION GRATING:     Ayush (992501030006)\n");
    printf("4. MALUS LAW & SPECIFIC ROTATION: Parth (992501030016)\n");
    printf("\nReferences:\n");
    printf("1. \"Optics\" by Ajoy Ghatak.\n");
    printf("2. Laboratory Manuals for Engineering Physics.\n");
    printf("3. Online resources on Newton's Rings, Diffraction,\n");
    printf("   and Wave Optics.\n");
    printf("=======================================================\n");
}
