using System;
using System.Numerics;

class Program
{
    static void Main()
    {
        // Define the parameters of the system
        int numPoints = 1000;                   // number of points to simulate
        double positionRange = 10.0;            // range of positions to simulate
        double mass = 1.0;                      // mass of the particle
        double potentialEnergy = 0.0;           // constant potential energy
        double dt = 0.001;                      // time step size
        double hbar = 1.0;                      // Planck's constant divided by 2*pi

        // Initialize the wavefunction
        Complex[] wavefunction = new Complex[numPoints];
        for (int i = 0; i < numPoints; i++)
        {
            double x = (i - numPoints / 2.0) * positionRange / numPoints;
            wavefunction[i] = Complex.Exp(-x * x / 2.0) / Math.Sqrt(Math.Sqrt(Math.PI));
        }

        // Evolve the wavefunction in time
        for (int step = 0; step < 10000; step++)
        {
            // Calculate the kinetic energy
            double kineticEnergy = 0.0;
            double k = 2.0 * Math.PI / (numPoints * positionRange);
            for (int i = 0; i < numPoints; i++)
            {
                double p = (i - numPoints / 2.0) * k;
                kineticEnergy += p * p * hbar * hbar / (2.0 * mass * numPoints * numPoints);
            }

            // Apply the time evolution operator
            Complex[] newWavefunction = new Complex[numPoints];
            for (int i = 0; i < numPoints; i++)
            {
                double x = (i - numPoints / 2.0) * positionRange / numPoints;
                Complex exponent = Complex.Exp(-Complex.ImaginaryOne * (kineticEnergy + potentialEnergy) * dt / hbar);
                newWavefunction[i] = wavefunction[i] * exponent;
            }
            wavefunction = newWavefunction;

            // Normalize the wavefunction
            double norm = 0.0;
            for (int i = 0; i < numPoints; i++)
            {
                norm += wavefunction[i].Magnitude * wavefunction[i].Magnitude;
            }
            for (int i = 0; i < numPoints; i++)
            {
                wavefunction[i] /= Math.Sqrt(norm);
            }

            // Print out the wavefunction every 100 steps
            if (step % 100 == 0)
            {
                Console.WriteLine($"Step {step}: {wavefunction[500].Real}");
            }
        }
    }
}
