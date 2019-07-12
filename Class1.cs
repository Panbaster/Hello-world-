using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace WindowsFormsApp3
{
    [Serializable()]
    class Network
    {
        double[][][,] network;
        double[][][,] neurons;
        public double alpha, LearnRate;
        readonly int EndNetworkLayers, NeuronX, NeuronY, NeuronsLayers;

        public Network(int InputsX, int InputsY, int NeuronX1, int NeuronY1, int NeuronsLayers1, int OutputCount, int EndNetworkLayers1 = 0, float alpha1 = 1, double LearnRate1 = 1)
        {
            EndNetworkLayers = EndNetworkLayers1;
            alpha = alpha1;
            NeuronX = NeuronX1;
            NeuronY = NeuronY1;
            LearnRate = LearnRate1;
            NeuronsLayers = NeuronsLayers1;

            double x, y;
            x = InputsX;
            y = InputsY;
            int k = 1;
            //Количество Слоев
            while (x >= NeuronX || y >= NeuronY)
            {
                k++;
                x = Math.Ceiling((x + 1 - NeuronX) / 2);
                y = Math.Ceiling((y + 1 - NeuronY) / 2);
            }
            network = new double[k + EndNetworkLayers + 1][][,];

            //Глубина слоев
            for (int i = 0; i < k; i++)
            {
                double j = 3 * Math.Pow(NeuronsLayers, i);
                network[i] = new double[Convert.ToInt32(j)][,];
            }


            x = InputsX;
            y = InputsY;
            k = 0;
            //Размерность матрицы слоя
            for (int i = 0; i < network[k].Length; i++)
            {
                network[k][i] = new double[Convert.ToInt32(x), Convert.ToInt32(y)];
            }
            x = x + 1 - NeuronX;
            y = y + 1 - NeuronY;
            k++;
            while (x > NeuronX || y > NeuronY)
            {

                x = Math.Ceiling(x / 2);
                y = Math.Ceiling(y / 2);
                for (int i = 0; i < network[k].Length; i++)
                {
                    network[k][i] = new double[Convert.ToInt32(x), Convert.ToInt32(y)];
                }
                k++;
                x = x + 1 - NeuronX;
                y = y + 1 - NeuronY;
            }
            for (int i = 0; i < network[network.Length - 2 - EndNetworkLayers].Length; i++)
            {
                network[network.Length - 2 - EndNetworkLayers][i] = new double[1, 1];
            }


            //Разметка Forward Propagation сети
            for (int i = 0; i < EndNetworkLayers; i++)
            {
                network[network.Length - 1 - EndNetworkLayers + i] = new double[network[network.Length - 2 - EndNetworkLayers + i].Length][,];
                for (int j = 0; j < network[network.Length - 2 - EndNetworkLayers + i].Length; j++)
                    network[network.Length - 1 - EndNetworkLayers + i][j] = new double[1, 1];
            }

            //Разметка вывода
            network[network.Length - 1] = new double[OutputCount][,];
            for (int i = 0; i < OutputCount; i++)
                network[network.Length - 1][i] = new double[1, 1];


            //Иницыацыя ядер
            neurons = new double[network.Length - 1][][,];
            Random rad = new Random();
            for (int i = 0; i < neurons.Length - 1 - EndNetworkLayers; i++)
            {
                neurons[i] = new double[NeuronsLayers][,];
                for (int j = 0; j < NeuronsLayers; j++)
                {
                    neurons[i][j] = new double[NeuronX, NeuronY];
                    for (int n = 0; n < NeuronX; n++)
                        for (int m = 0; m < NeuronY; m++)
                            neurons[i][j][n, m] = (rad.NextDouble() * 2) - 1;
                }
            }
            //Forward propagation's weights
            for (int i = neurons.Length - 1 - EndNetworkLayers; i < neurons.Length - 1; i++)
            {
                neurons[i] = new double[Convert.ToInt32(Math.Pow(network[network.Length - 2].Length, 2))][,];
                for (int j = 0; j < neurons[i].Length; j++)
                {
                    neurons[i][j] = new double[1, 1];
                    neurons[i][j][0, 0] = (rad.NextDouble() * 2) - 1;
                }
            }
            neurons[neurons.Length - 1] = new double[OutputCount * network[network.Length - 2].Length][,];
            for (int i = 0; i < neurons[neurons.Length - 1].Length; i++)
            {
                neurons[neurons.Length - 1][i] = new double[1, 1];
                neurons[neurons.Length - 1][i][0, 0] = Convert.ToDouble((rad.NextDouble() * 2) - 1);
            }

        }

        private double[,] ArrMultiplay(double[,] neuron, double[,] arr)
        {
            double[,] Out = new double[arr.GetLength(0) - neuron.GetLength(0) + 1, arr.GetLength(1) - neuron.GetLength(1) + 1];
            double[,] Out2 = new double[Convert.ToInt32(Math.Ceiling(Out.GetLength(0) / 2D)), Convert.ToInt32(Math.Ceiling(Out.GetLength(1) / 2D))];

            for (int i = 0; i < Out2.GetLength(1); i++)
                for (int j = 0; j < Out2.GetLength(0); j++)
                    Out2[j, i] = -1D;

            for (int i = 0; i < Out.GetLength(1); i++)
                for (int j = 0; j < Out.GetLength(0); j++)
                {
                    for (int n = 0; n < neuron.GetLength(1); n++)
                        for (int m = 0; m < neuron.GetLength(0); m++)
                        {
                            Out[j, i] += neuron[m, n] * arr[j + m, i + n];
                        }
                    Out[j, i] = 2D / (1D + Math.Pow(Math.E, -Out[j, i] * alpha)) - 1D;       //функцыя активацыы
                }


            for (int i = 0; i < Out.GetLength(1); i++)      //pulling
                for (int j = 0; j < Out.GetLength(0); j++)
                    if (Out2[Convert.ToInt32(Math.Floor(j / 2D)), Convert.ToInt32(Math.Floor(i / 2D))] < Out[j, i])
                        Out2[Convert.ToInt32(Math.Floor(j / 2D)), Convert.ToInt32(Math.Floor(i / 2D))] = Out[j, i];

            if (Out2.GetLength(0) < neuron.GetLength(0) & Out2.GetLength(1) < neuron.GetLength(1)) //flooring for last layer
            {
                double a = -1;
                for (int i = 0; i < Out2.GetLength(1); i++)
                    for (int j = 0; j < Out2.GetLength(0); j++)
                        if (a < Out2[j, i])
                            a = Out2[j, i];
                Out2 = new double[1, 1];
                Out2[0, 0] = a;
            }

            return Out2;
        }

        private double[][,] NetworkResolve(double[][,] neuron, double[][,] arr)
        {
            double[][,] Out = new double[neuron.Length / arr.Length][,];
            for (int i = 0; i < Out.Length; i++)
            {
                Out[i] = new double[1, 1];
                for (int j = 0; j < arr.Length; j++)
                    Out[i][0, 0] += arr[j][0, 0] * neuron[(arr.Length * i) + j][0, 0];

                Out[i][0, 0] = 2 / (1 + Math.Pow(Math.E, -Out[i][0, 0] * alpha)) - 1;      //функцыя активацыы
            }
            return Out;
        }

        public double[] Compute(double[][,] Input)
        {
            if (Verify(Input))
                throw new System.ArgumentException("Input size (" + Input[0].GetLength(0) + " " + Input[0].GetLength(1)+")wrong. suposed to be " + network[0].GetLength(0) + " " + network[0].GetLength(0));

                network[0] = Input;
            double[] Out = new double[network[network.Length - 1].Length];
           

            for (int i = 1; i < network.Length - 1 - EndNetworkLayers; i++)         //расчет сверточной сети
                for (int j = 0; j < network[i - 1].Length; j++)
                    for (int n = 0; n < NeuronsLayers; n++)
                        network[i][j * NeuronsLayers + n] = ArrMultiplay(neurons[i - 1][n], network[i - 1][j]);
            

            for (int i = network.Length - EndNetworkLayers - 2; i < network.Length - 1; i++) //расчет полносвязной сети
                network[i + 1] = NetworkResolve(neurons[i], network[i]);
            

            for (int i = 0; i < network[network.Length - 1].Length; i++)
                Out[i] = network[network.Length - 1][i][0, 0];

            return Out;

        }

        private double[,] ConvErrorResolve(double[,] neuron, double[,] arr, int X, int Y)
        {
            double[,] arr2 = new double[X - neuron.GetLength(0) + 1, Y - neuron.GetLength(1) + 1];
            for (int i = 0; i < arr2.GetLength(1); i++)
                for (int j = 0; j < arr2.GetLength(0); j++)
                    arr2[i, j] = arr[Convert.ToInt32(Math.Floor(i / 2D)), Convert.ToInt32(Math.Floor(j / 2D))];

            double[,] Out = new double[X, Y];
            for (int i = 0; i < Out.GetLength(1); i++)
                for (int j = 0; j < Out.GetLength(0); j++)
                {
                    if (j >= 60)
                    { }
                    for (int n = Math.Max(0, i - Out.GetLength(1) + neuron.GetLength(1)); n < neuron.GetLength(1) & n <= i; n++)
                        for (int m = Math.Max(0, j - Out.GetLength(0) + neuron.GetLength(0)); m < neuron.GetLength(0) & m <= j; m++)
                        {
                            Out[j, i] += neuron[m, n] * arr2[j - m, i - n];
                        }
                }
            return Out;
        }

        private bool Verify(double[][,] Input)
        {
            bool Out = false;
            for (int i = 0; i < Input.Length; i++)
                if (Input[i].GetLength(0) != network[0][i].GetLength(0) & Input[i].GetLength(1) != network[0][i].GetLength(1))
                    Out = true;
            return Out;
        }

        public double[] Teach(double[][,] Input, double[] Ansver)
        {
            if (Verify(Input) || Ansver.Length != network[network.Length-1].Length)
                throw new System.ArgumentException("Input size (" + Input[0].GetLength(0) + " " + Input[0].GetLength(1)+")wrong. suposed to be " + network[0].GetLength(0) + " " + network[0].GetLength(0));

            double[] ErrorOut = new double[network[network.Length - 1].Length];
            double[][][,] error = new double[network.Length][][,];      //разметка памяти под масив ошыбок
            for (int i = 0; i < network.Length; i++)
            {
                error[i] = new double[network[i].Length][,];
                for (int j = 0; j < network[i].Length; j++)
                    error[i][j] = new double[network[i][j].GetLength(0), network[i][j].GetLength(1)];
            }
            Compute(Input);


            for (int i = 0; i < network[network.Length - 1].Length; i++)      //расчет ошыбки вывода сети
            {
                error[network.Length - 1][i][0, 0] = Ansver[i] - network[network.Length - 1][i][0, 0];
                ErrorOut[i] = Ansver[i] - network[network.Length - 1][i][0, 0];
            }

            for (int i = network.Length - 2; i >= network.Length - 2 - EndNetworkLayers; i--)       //расчет ошыбки полносвязной сети
                for (int j = 0; j < error[i + 1].Length; j++)
                    for (int n = 0; n < error[i].Length; n++)
                        error[i][n][0, 0] += error[i + 1][j][0, 0] * neurons[i][n + (j * error[i].Length)][0, 0];

            double[,] buff = new double[Convert.ToInt32(Math.Ceiling((network[network.Length - EndNetworkLayers - 3][0].GetLength(0) + 1 - NeuronX) / 2D)),
                Convert.ToInt32(Math.Ceiling((network[network.Length - EndNetworkLayers - 3][0].GetLength(1) + 1 - NeuronY) / 2D))];
            double[][,] buff2 = new double[NeuronsLayers][,];
            for (int i = 0; i < network[network.Length - 3 - EndNetworkLayers].Length; i++)    //заполнение ошыбок первого слоя
            {
                for (int j = 0; j < NeuronsLayers; j++)
                {
                    for (int n = 0; n < buff.GetLength(1); n++)     //генерацыя стартовой матрицы
                        for (int m = 0; m < buff.GetLength(0); m++)
                            buff[m, n] = error[network.Length - 2 - EndNetworkLayers][(i * NeuronsLayers) + j][0, 0];

                    buff2[j] = ConvErrorResolve(neurons[network.Length - 3 - EndNetworkLayers][j], buff,    //буферизацыя матриц для сложения
                        error[network.Length - 3 - EndNetworkLayers][i].GetLength(0), error[network.Length - 3 - EndNetworkLayers][i].GetLength(1));

                }
                for (int j = 0; j < NeuronsLayers; j++)
                    for (int n = 0; n < buff2[j].GetLength(1); n++)
                        for (int m = 0; m < buff2[j].GetLength(0); m++)
                            error[network.Length - 3 - EndNetworkLayers][i][m, n] += buff2[j][m, n];
            }
            for (int k = network.Length - 4 - EndNetworkLayers; k >= 0; k--)
                for (int i = 0; i < network[k].Length; i++) //заполнение ошыбок сверточной сети
                {
                    for (int j = 0; j < NeuronsLayers; j++)
                        buff2[j] = ConvErrorResolve(neurons[k][j], network[k + 1][j + (i * NeuronsLayers)],    //буферизацыя матриц для сложения
                            error[k][i].GetLength(0), error[k][i].GetLength(1));

                    for (int j = 0; j < NeuronsLayers; j++)
                        for (int n = 0; n < buff2[j].GetLength(1); n++)
                            for (int m = 0; m < buff2[j].GetLength(0); m++)
                                error[k][i][m, n] += buff2[j][m, n];
                }

            //------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            int a = 0;
            for (int i = network.Length - 1; i > network.Length - 2 - EndNetworkLayers; i--)   //поправка весов полносвязной сети
                for (int j = 0; j < neurons[i - 1].Length; j++)
                {
                    Math.DivRem(j, network[i - 1].Length, out a);
                    neurons[i - 1][j][0, 0] += alpha *
                        Math.Abs(network[i - 1][a][0, 0]) * //выход преведушего нейрона
                        error[i][Convert.ToInt32(Math.Floor((double)j / network[i - 1].Length))][0, 0] *    //ошыбка следуюшего нейрона
                        network[i][Convert.ToInt32(Math.Floor((double)j / network[i - 1].Length))][0, 0] * (1 - network[i][Convert.ToInt32(Math.Floor((double)j / network[i - 1].Length))][0, 0]);  //производная следуйшенго нейрона
                }

            //поправка весов сверточной сети
            buff = new double[error[error.Length - 3 - EndNetworkLayers][0].GetLength(0) - NeuronX + 1, error[error.Length - 3 - EndNetworkLayers][0].GetLength(1) - NeuronY + 1];
            double[,] buff3 = new double[buff.GetLength(0), buff.GetLength(1)];
            for (int n = 0; n < network[network.Length - 2 - EndNetworkLayers].Length; n++)     //обработка первого слоя ядер
            {
                for (int i = 0; i < buff3.GetLength(1); i++)
                    for (int j = 0; j < buff3.GetLength(0); j++)
                        buff3[j, i] = network[network.Length - 2 - EndNetworkLayers][n][0, 0];
                for (int i = 0; i < buff.GetLength(1); i++)
                    for (int j = 0; j < buff.GetLength(0); j++)
                        buff[j, i] = error[network.Length - 2 - EndNetworkLayers][n][0, 0];

                int b = 0;
                b = Math.DivRem( n , neurons[neurons.Length - 2 - EndNetworkLayers].Length, out a);
                neurons[neurons.Length - 2 - EndNetworkLayers][a] = ConvErrorResolve(neurons[neurons.Length - 2 - EndNetworkLayers][a], network[network.Length - 3 - EndNetworkLayers][b], buff, buff3);
            }
            for (int i = network.Length - 4 - EndNetworkLayers; i >= 0; i--) //Обработка остальных слоев
                for (int n = 0; n < network[i+1].Length; n++)
                {
                    int b = Math.DivRem(n, neurons[neurons.Length - 2 - EndNetworkLayers].Length, out a);
                    neurons[i][a] = ConvErrorResolve(neurons[i][a], network[i][b], error[i+1][n], network[i+1][n]);
                }

            return ErrorOut;
        }

        double [,] ConvErrorResolve(double[,] neuron, double[,] arr,double[,] error,double[,] derivative)
        {
            double[,] Out = neuron;
            for (int i=0;i<derivative.GetLength(1);i++)
                for (int j=0;j<derivative.GetLength(0);j++)
                    for (int n=0;n<neuron.GetLength(1);n++)
                        for (int m=0;m<neuron.GetLength(0);m++)
                            neuron[m,n] += LearnRate*
                                error[j,i]*
                                Math.Abs(arr[m,n])*
                                (2*Math.Pow(Math.E, -derivative[j, i]))/Math.Pow(2,1+Math.Pow(Math.E,-derivative[j, i])); //Производная функцыы активацыы
            return Out;
        }

    }
}