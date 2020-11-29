using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Com_Methods
{
    class CG
    {
        // макс. итераций
        public int Max_Iter { set; get; }

        // точность невязки
        public double Eps { set; get; }

        // текущая итерация
        public int Iter { set; get; }

        // предобусловливатель
        public IPredconditioner Pr { set; get; }

        public CG(int MAX, double EPS)
        {
            Max_Iter = MAX;
            Eps = EPS;
            Iter = 0;
        }

        public Vector Start_Solver(CSlR_Matrix A, Vector F, int PREC)
        {
            int n = A.N;
            // решение
            Vector RES = new Vector(n);
            for (int i = 0; i < n; i++)
                RES.Elem[i] = 0.0;
            //реализация предобусловливателя
            switch (PREC)
            {
                case 1: { Pr = new Diag_Predconditioner(A); break; }
                case 2: { Pr = new ILU_Decomposition(A); break; }
                //default: { Pr = new E_Predconditioner(); break; }
            }
            var r = new Vector(n);
            var p = new Vector(n);
            var vec = new Vector(n);

            // параметры 
            double alpha, betta, sc1, sc2;

            // флаг для остановки
            int Flag = 0;

            // ||r||
            double nev_r = 0;

            // начало CG 


            A.Mult_MV(RES, vec);
            for (int i = 0; i < n; i++)
            {
                r.Elem[i] = F.Elem[i] - vec.Elem[i];
            }
            Pr.Start_Predconditioner(r, p); // в r предобусловленная система-невязка

            while (Flag != 1 && Iter < Max_Iter)
            {
                sc1 = p * r;
                A.Mult_MV(p, vec);
                sc2 = vec * p;
                alpha = sc1 / sc2;
                for (int i = 0; i < n; i++)
                {
                    RES.Elem[i] += alpha * p.Elem[i];
                    r.Elem[i] -= alpha * vec.Elem[i];
                }
                Pr.Start_Tr_Predconditioner(r, vec);

                sc2 = vec * p;
                betta = sc2 / sc1;
                nev_r = r.Norm();
                if (nev_r < Eps)
                    Flag = 1;
                if (Flag == 0)
                {
                    for (int i = 0; i < n; i++)
                    {
                        p.Elem[i] = vec.Elem[i] + betta * p.Elem[i];
                    }
                }
                Console.WriteLine("{0,-20} {1,-20}", Iter, nev_r.ToString("E"));
                Iter++;
            }
            return RES;
        }
    }
}
