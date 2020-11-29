using System;
using System.Diagnostics;
using System.Threading;


namespace Com_Methods
{
    class CONST
    {
        //точность решения
        public static double EPS = 1e-16;
    }

    class Tools
    {
        //замер времени
        public static string Measurement_Time(Thread thread)
        {
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();
            thread.Start();
            while (thread.IsAlive) ;
            stopWatch.Stop();
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
            return ("RunTime: " + elapsedTime);
        }

        //относительная погрешность
        public static double Relative_Error(Vector X, Vector x)
        {
            double s = 0.0;
            for (int i = 0; i < X.N; i++)
            {
                s += Math.Pow(X.Elem[i] - x.Elem[i], 2);
            }
            return Math.Sqrt(s) / x.Norm();
        }

        //относительная невязка
        public static double Relative_Discrepancy (Matrix A, Vector X, Vector F)
        {
            var x = A * X;
            for (int i = 0; i < X.N; i++) x.Elem[i] -= F.Elem[i];
            return x.Norm() / F.Norm();
        }
    }

    class Program
    {
        static void Main()
        {
            try
            {
                //методы на подпространствах Крылова: CSlR-формат матрицы
                var T3 = new Thread(() =>
                {
                    //разреженная матрица
                    var A = new CSlR_Matrix("Data\\Sparse_Format\\Systems1\\NonSPD\\");
                    //заполнение вектора истинного решения и правой части
                    var F = new Vector(A.N);
                    var X_True = new Vector(A.N);
                    int pre = 1;
                    //заполнение вектора истинного решения и правой части
                    for (int i = 0; i < A.N; i++) { X_True.Elem[i] = 1.0; }
                    A.Mult_MV(X_True, F);
                    var Solver = new BCG(30000, 1e-12);
                    var X = Solver.Start_Solver(A, F, 2);
                    Console.WriteLine(Tools.Relative_Error(X, X_True));
                });

                //время решения
                Console.WriteLine(Tools.Measurement_Time(T3));
                
            }
            catch (Exception E)
            {
                Console.WriteLine("\n*** Error! ***");
                Console.WriteLine("Method:  {0}",   E.TargetSite);
                Console.WriteLine("Message: {0}\n", E.Message);
            }
        }
    }
}