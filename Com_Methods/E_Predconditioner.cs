using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Com_Methods
{
    public interface IPredconditioner
    {
        //реализация
        void Start_Predconditioner(Vector X, Vector RES);

        //Реализация для Aт
        void Start_Tr_Predconditioner(Vector X, Vector RES);
    }

    class E_Predconditioner
    {
        public E_Predconditioner() { }
        public void Start_Predconditioner(Vector X, Vector RES)
        {
            for (int i = 0; i < X.N; i++)
                RES.Elem[i] = X.Elem[i];
        }

        public void Start_Tr_Predconditioner(Vector X, Vector RES)
        {
            for (int i = 0; i < X.N; i++)
                RES.Elem[i] = X.Elem[i];
        }
    }
    class Diag_Predconditioner : IPredconditioner
    {
        Vector Diag { set; get; }
        public Diag_Predconditioner(CSlR_Matrix A)
        {
            Diag = new Vector(A.N);
            for (int i = 0; i < A.N; i++)
                Diag.Elem[i] = A.di[i];
        }

        public void Start_Predconditioner(Vector X, Vector RES)
        {
            for (int i = 0; i < X.N; i++)
                RES.Elem[i] = X.Elem[i] / Diag.Elem[i];
        }

        public void Start_Tr_Predconditioner(Vector X, Vector RES)
        {
            for (int i = 0; i < X.N; i++)
                RES.Elem[i] = X.Elem[i] / Diag.Elem[i];
        }
    }

    class ILU_Decomposition : IPredconditioner
    {
        //хранилище для LU-разложения
        CSlR_Matrix LU { set; get; }
        public ILU_Decomposition(CSlR_Matrix A)
        {
            LU =A.Create_ILU_Decomposition();
        }

        public void Start_Predconditioner(Vector X, Vector RES)
        {
            //подкаст: интерфейс приводится к типу
            (LU as CSlR_Matrix).Slau_L(RES, X);
            (LU as CSlR_Matrix).Slau_U(RES, RES);
        }
        public void Start_Tr_Predconditioner(Vector X, Vector RES)
        {
            //подкаст: интерфейс приводится к типу
            (LU as CSlR_Matrix).Slau_Ut(RES, X);
            (LU as CSlR_Matrix).Slau_Lt(RES, RES);
        }

    }
}
