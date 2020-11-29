using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Com_Methods
{
    class BCG
    {
        // макс. итераций
        public int Max_Iter { set; get; }

        // точность невязки
        public double Eps { set; get; }

        // текущая итерация
        public int Iter { set; get; }

        // предобусловливатель
        public IPredconditioner Pr { set; get; }

        public BCG(int MAX, double EPS)
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
               // default: { Pr = new E_Predconditioner(); break; }
            }
            var r = new Vector(n);
            var p = new Vector(n);
            var r_ = new Vector(n);
            var p_ = new Vector(n);
            var vec = new Vector(n);
            var vec_ = new Vector(n);
            var vec_help = new Vector(n);

            // параметры 
            double alpha, betta, sc1, sc2;

            // флаг для остановки
            int Flag = 0;

            // ||r||
            double nev_r = 0;

            // начало Bi_Conjugate_Gradient_Method 


            A.Mult_MV(RES, vec_help);
            for (int i = 0; i < n; i++)
            {
                vec_help.Elem[i] = F.Elem[i] - vec_help.Elem[i];
            }
            Pr.Start_Predconditioner(vec_help, r); // в r предобусловленная система-невязка
            p.Copy(r);
            r_.Copy(r);
            p_.Copy(r_);
            while (Flag != 1 && Iter < Max_Iter)
            {
                sc1 = r * r_;
                A.Mult_MV(p, vec_help);
                Pr.Start_Predconditioner(vec_help, vec);
                sc2 = vec * p_;
                alpha = sc1 / sc2;
                for (int i = 0; i < n; i++)
                {
                    RES.Elem[i] += alpha * p.Elem[i];
                    r.Elem[i] -= alpha * vec.Elem[i];
                }
                Pr.Start_Tr_Predconditioner(p_, vec_help);
                A.Mult_MtV(vec_help, vec_);
                for (int i = 0; i < n; i++)
                {
                    r_.Elem[i] -= alpha * vec_.Elem[i];
                }
                sc2 = r * r_;
                betta = sc2 / sc1;
                nev_r = r.Norm();
                if (nev_r < Eps)
                    Flag = 1;
                if (Flag == 0)
                {
                    for (int i = 0; i < n; i++)
                    {
                        p.Elem[i] = r.Elem[i] + betta * p.Elem[i];
                        p_.Elem[i] = r_.Elem[i] + betta * p_.Elem[i];
                    }
                }
                //Console.WriteLine("{0,-20} {1,-20}", Iter, nev_r.ToString("E"));
                Console.WriteLine(nev_r.ToString("E"));
                Iter++;
            }
            return RES;
        }

    }
}
