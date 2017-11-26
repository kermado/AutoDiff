using System;
using System.Diagnostics;

namespace AutoDiff
{
    public struct Jet
    {
        private double r;
        private double[] v;

        public double Real { get { return this.r; } }
        public double[] Infinitesimals { get { return this.v; } }

        public Jet(double value, int index, int dimension)
        {
            this.r = value;
            this.v = new double[dimension];
            this.v[index] = 1.0;
        }

        private Jet(double r, double[] v)
        {
            this.r = r;
            this.v = v;
        }

        #region Array_Operations

        private static double[] Negative(double[] a)
        {
            return Negative(a, a.Length);
        }

        private static double[] Negative(double[] a, int dimension)
        {
            double[] result = new double[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                result[i] = -a[i];
            }

            return result;
        }

        private static double[] Add(double[] a, double[] b)
        {
            Debug.Assert(a.Length == b.Length);
            return Add(a, b, a.Length);
        }

        private static double[] Add(double[] a, double[] b, int dimension)
        {
            double[] result = new double[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                result[i] = a[i] + b[i];
            }

            return result;
        }

        private static double[] Add(double[] a, double s)
        {
            return Add(a, s, a.Length);
        }

        private static double[] Add(double[] a, double s, int dimension)
        {
            double[] result = new double[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                result[i] = a[i] + s;
            }

            return result;
        }

        private static double[] Subtract(double[] a, double[] b)
        {
            Debug.Assert(a.Length == b.Length);
            return Subtract(a, b, a.Length);
        }

        private static double[] Subtract(double[] a, double[] b, int dimension)
        {
            double[] result = new double[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                result[i] = a[i] - b[i];
            }

            return result;
        }

        private static double[] Subtract(double[] a, double s)
        {
            return Subtract(a, s, a.Length);
        }

        private static double[] Subtract(double[] a, double s, int dimension)
        {
            double[] result = new double[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                result[i] = a[i] - s;
            }

            return result;
        }

        private static double[] Multiply(double[] a, double[] b)
        {
            Debug.Assert(a.Length == b.Length);
            return Multiply(a, b, a.Length);
        }

        private static double[] Multiply(double[] a, double[] b, int dimension)
        {
            double[] result = new double[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                result[i] = a[i] * b[i];
            }

            return result;
        }

        private static double[] Multiply(double[] a, double s)
        {
            return Multiply(a, s, a.Length);
        }

        private static double[] Multiply(double[] a, double s, int dimension)
        {
            double[] result = new double[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                result[i] = a[i] * s;
            }

            return result;
        }

        private static double[] Divide(double[] a, double[] b)
        {
            Debug.Assert(a.Length == b.Length);
            return Divide(a, b, a.Length);
        }

        private static double[] Divide(double[] a, double[] b, int dimension)
        {
            double[] result = new double[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                result[i] = a[i] / b[i];
            }

            return result;
        }

        private static double[] Divide(double[] a, double s)
        {
            return Divide(a, s, a.Length);
        }

        private static double[] Divide(double[] a, double s, int dimension)
        {
            double[] result = new double[dimension];
            for (int i = 0; i < dimension; ++i)
            {
                result[i] = a[i] / s;
            }

            return result;
        }

        #endregion

        #region Binary_Addition

        public static Jet operator +(Jet a, Jet b)
        {
            return new Jet(a.r + b.r, Add(a.v, b.v));
        }

        public static Jet operator +(Jet a, double s)
        {
            return new Jet(a.r + s, a.v);
        }

        public static Jet operator +(double s, Jet a)
        {
            return new Jet(s + a.r, a.v);
        }

        #endregion

        #region Binary_Subtraction

        public static Jet operator -(Jet a, Jet b)
        {
            return new Jet(a.r - b.r, Subtract(a.v, b.v));
        }

        public static Jet operator -(Jet a, double s)
        {
            return new Jet(a.r - s, a.v);
        }

        public static Jet operator -(double s, Jet a)
        {
            return new Jet(s - a.r, a.v);
        }

        #endregion

        #region Binary_Multiplication

        public static Jet operator *(Jet a, Jet b)
        {
            Debug.Assert(a.v.Length == b.v.Length);
            int n = a.v.Length;
            double[] v = new double[n];
            for (int i = 0; i < n; ++i)
            {
                v[i] = a.v[i] * b.r + b.v[i] * a.r;
            }

            return new Jet(a.r * b.r, v);
        }

        public static Jet operator *(Jet a, double s)
        {
            return new Jet(a.r * s, Multiply(a.v, s));
        }

        public static Jet operator *(double s, Jet a)
        {
            return new Jet(a.r * s, Multiply(a.v, s));
        }

        #endregion

        #region Binary_Division

        public static Jet operator /(Jet a, Jet b)
        {
            double t1 = 1.0 / b.r;
            double t2 = a.r * t1;

            Debug.Assert(a.v.Length == b.v.Length);
            int n = a.v.Length;
            double[] v = new double[n];
            for (int i = 0; i < n; ++i)
            {
                v[i] = (a.v[i] - t2 * b.v[i]) * t1;
            }

            return new Jet(a.r * t1, v);
        }

        public static Jet operator /(Jet a, double s)
        {
            double t = 1.0 / s;
            return new Jet(a.r * t, Multiply(a.v, t));
        }

        public static Jet operator /(double s, Jet a)
        {
            double t = -s / (a.r * a.r);
            return new Jet(s / a.r, Multiply(a.v, t));
        }

        #endregion

        #region Unary

        public static Jet operator +(Jet a)
        {
            return a;
        }

        public static Jet operator -(Jet a)
        {
            return new Jet(-a.r, Negative(a.v));
        }

        #endregion

        #region Comparison

        public static bool operator <(Jet a, Jet b)
        {
            return a.r < b.r;
        }

        public static bool operator <=(Jet a, Jet b)
        {
            return a.r <= b.r;
        }

        public static bool operator >(Jet a, Jet b)
        {
            return a.r > b.r;
        }

        public static bool operator >=(Jet a, Jet b)
        {
            return a.r >= b.r;
        }

        public static bool operator <(Jet a, double s)
        {
            return a.r < s;
        }

        public static bool operator <=(Jet a, double s)
        {
            return a.r <= s;
        }

        public static bool operator >(Jet a, double s)
        {
            return a.r > s;
        }

        public static bool operator >=(Jet a, double s)
        {
            return a.r >= s;
        }

        public static bool operator <(double s, Jet b)
        {
            return s < b.r;
        }

        public static bool operator <=(double s, Jet b)
        {
            return s <= b.r;
        }

        public static bool operator >(double s, Jet b)
        {
            return s > b.r;
        }

        public static bool operator >=(double s, Jet b)
        {
            return s >= b.r;
        }

        #endregion

        #region Functions

        public static double Abs(double s)
        {
            return Math.Abs(s);
        }

        public static double Log(double s)
        {
            return Math.Log(s);
        }

        public static double Exp(double s)
        {
            return Math.Exp(s);
        }

        public static double Sqr(double s)
        {
            return s * s;
        }

        public static double Sqrt(double s)
        {
            return Math.Sqrt(s);
        }

        public static double Sin(double s)
        {
            return Math.Sin(s);
        }

        public static double Cos(double s)
        {
            return Math.Cos(s);
        }

        public static double Tan(double s)
        {
            return Math.Tan(s);
        }

        public static double Asin(double s)
        {
            return Math.Asin(s);
        }

        public static double Acos(double s)
        {
            return Math.Acos(s);
        }

        public static double Atan(double s)
        {
            return Math.Atan(s);
        }

        public static Jet Abs(Jet a)
        {
            return (a.r < 0.0) ? -a : a;
        }

        public static Jet Log(Jet a)
        {
            double t = 1.0 / a.r;
            return new Jet(Log(a.r), Multiply(a.v, t));
        }

        public static Jet Exp(Jet a)
        {
            double t = Exp(a.r);
            return new Jet(t, Multiply(a.v, t));
        }

        public static Jet Sqr(Jet a)
        {
            return new Jet(Sqr(a.r), Multiply(a.v, (2.0 * a.r)));
        }

        public static Jet Sqrt(Jet a)
        {
            double t = Sqrt(a.r);
            return new Jet(t, Multiply(a.v, 1.0 / (2.0 * t)));
        }

        public static Jet Sin(Jet a)
        {
            return new Jet(Sin(a.r), Multiply(a.v, Cos(a.r)));
        }

        public static Jet Cos(Jet a)
        {
            return new Jet(Cos(a.r), Multiply(a.v, -Sin(a.r)));
        }

        public static Jet Tan(Jet a)
        {
            double t = Tan(a.r);
            return new Jet(t, Multiply(a.v, 1.0 + t + t));
        }

        public static Jet Asin(Jet a)
        {
            return new Jet(Asin(a.r), Multiply(a.v, 1.0 / Sqrt(1.0 - a.r * a.r)));
        }

        public static Jet Acos(Jet a)
        {
            return new Jet(Acos(a.r),  Multiply(a.v, -1.0 / Sqrt(1.0 - a.r * a.r)));
        }

        public static Jet Atan(Jet a)
        {
            return new Jet(Atan(a.r), Multiply(a.v, 1.0 / (1.0 + a.r * a.r)));
        }

        #endregion

        #region Classification

        public static bool IsInfinity(double s)
        {
            return Double.IsInfinity(s);
        }

        public static bool IsInfinity(Jet a)
        {
            if (IsInfinity(a.r)) { return true; }
            for (int i = 0; i < a.v.Length; ++i)
            {
                if (IsInfinity(a.v[i])) { return true; }
            }

            return false;
        }

        public static bool IsNaN(double s)
        {
            return Double.IsNaN(s);
        }

        public static bool IsNaN(Jet a)
        {
            if (IsNaN(a.r)) { return true; }
            for (int i = 0; i < a.v.Length; ++i)
            {
                if (IsNaN(a.v[i])) { return true; }
            }

            return false;
        }

        #endregion

        public override string ToString()
        {
            return String.Format("({0}, [{1}])", this.r, String.Join(", ", this.v));
        }
    }
}
