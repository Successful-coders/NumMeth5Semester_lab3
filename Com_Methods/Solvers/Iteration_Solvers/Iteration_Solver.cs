using System;

namespace Com_Methods
{
    public interface IIteration_Solver
    {
        //максимальное число итераций
        int Max_Iter { set; get; }
        //точность решения
        double Eps { set; get; }
        //текущая итерация 
        int Iter { set; get; }
    }
}
