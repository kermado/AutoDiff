using System;

namespace AutoDiff
{
    public delegate void DifferentiableFunction(Jet[] values, Jet[] parameters);

    public class AutoDiff
    {
        private int range;
        private int domain;
        private int index;
        private string[] parameters;

        public AutoDiff(int range, int domain)
        {
            this.range = range;
            this.domain = domain;
            this.index = 0;
            this.parameters = new string[domain];
        }

        public Jet Add(string parameter, double value = 0.0)
        {
            if (this.index >= domain) { throw new Exception("Too many parameters"); }
            int index = this.index++;
            this.parameters[index] = parameter;
            return new Jet(value, index, this.domain);
        }

        public void Evaluate(DifferentiableFunction function, Jet[] parameters, double[] values, double[][] jacobian)
        {
            if (this.index < this.domain) { throw new Exception("Insufficient parameters added"); }
            if (this.domain != parameters.Length) { throw new ArgumentException("Insufficient parameters provided"); }

            var f = new Jet[this.range];
            function(f, parameters);

            for (int i = 0; i < this.range; ++i)
            {
                values[i] = f[i].Real;

                for (int j = 0; j < this.domain; ++j)
                {
                    jacobian[i][j] = f[i].Infinitesimals[j];
                }
            }
        }
    }
}
