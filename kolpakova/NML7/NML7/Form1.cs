using System;
using System.Windows.Forms;
using System.Collections.Generic;

namespace NML5
{
    public partial class Form1 : Form
    {
        
        int M;
        static int N = 10;
        int K = N;
        double l = 1;
        double[,] U = new double[N + 1, N + 1];
        double[] dt = new double[N + 1];
        double[,] U2 = new double[N + 1, N + 1];
        double[] dt2 = new double[N + 1];
        double[,] U3 = new double[N + 1, N + 1];
        double[] dt3 = new double[N + 1];
        static int Zcount = 0;
        static int Pcount = 0;
        static int PUcount = 0;
        public Form1()
        {
            InitializeComponent();
            trackBar1.Maximum = K;
            trackBar1.Minimum = -K;
        }
        private void paint()
        {
            double h = l / N;
            double y = 0, x = 0;
            double u = M / Convert.ToDouble(K);
            chart1.Series[0].Points.Clear();
            chart1.ChartAreas[0].AxisY.LabelStyle.Format = "f2";
            chart1.ChartAreas[0].AxisX.LabelStyle.Format = "f2";
            chart1.ChartAreas[0].AxisX.Maximum = l;
            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisY.Maximum = l; 
            chart1.ChartAreas[0].AxisY.Minimum = 0;
            if (u >= 0)
            {
                for (double i = 0; i <= l; i += h)
                {
                    y = i;
                    x = Math.Sqrt(Math.Pow(y, 2) + u);
                    chart1.Series[0].Points.AddXY(x, y);
                }
            }
            else
            {
                for (double i = 0; i <= l; i += h)
                {
                    x = i;
                    y = Math.Sqrt(Math.Pow(x, 2) - u);
                    chart1.Series[0].Points.AddXY(x, y);
                }
            }
        }
        private void paintY()
        {
            double h = l / N;
            double u = 0, x = 0;
            double y = M / Convert.ToDouble(K);
            chart1.Series[0].Points.Clear();
            chart1.ChartAreas[0].AxisY.LabelStyle.Format = "f2";
            chart1.ChartAreas[0].AxisX.LabelStyle.Format = "f2";
            chart1.ChartAreas[0].AxisX.Maximum = l;
            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisY.Maximum = l;
            chart1.ChartAreas[0].AxisY.Minimum = 0;
            for (double i = 0; i <= l; i += h)
            {
                x = i;
                u = Math.Pow(x, 2) - Math.Pow(y, 2);
                chart1.Series[0].Points.AddXY(x, u);
            }
        }
        private void paintX()
        {
            double h = l / N;
            double y = 0, u = 0;
            double x = M / Convert.ToDouble(K);
            chart1.Series[0].Points.Clear();
            chart1.ChartAreas[0].AxisY.LabelStyle.Format = "f2";
            chart1.ChartAreas[0].AxisX.LabelStyle.Format = "f2";
            chart1.ChartAreas[0].AxisX.Maximum = l;
            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisY.Maximum = l;
            chart1.ChartAreas[0].AxisY.Minimum = 0;
            for (double i = 0; i <= l; i += h)
            {
                y = i;
                u = Math.Pow(x, 2) - Math.Pow(y, 2);
                chart1.Series[0].Points.AddXY(y, u);
            }
        }
        private void CRS()
        {
            double h = l / N;
            double[,] u = new double[K + 1, N + 1];
            double[,] matr = new double[(N + 1) * (N + 1), (N + 1) * (N + 1)];
            double[] d = new double[(N + 1) * (N + 1)];
            for (int j = 0; j <= N - 1; j++)
            {
                matr[j, j] = -1; matr[j, N + 1 + j] = 1;
            }
            matr[N, N] = 1; d[N] = - 1;
            matr[(N + 1) * (N + 1) - 1, (N + 1) * (N + 1) - 1] = 1; d[(N + 1) * (N + 1) - 1] = 0;
            for (int j = (N+1)*(N+1)-2; j > (N + 1) * (N + 1) - N - 2; j--)
            {
                matr[j, j] = 1;
                d[j] = 1 - Math.Pow(h * (j - (N + 1) * (N + 1) + N + 1), 2);
            }
            for (int i = 1; i <= N - 1; i++)
            {
                matr[(N + 1) * i, (N + 1) * i] = -1;
                matr[(N + 1) * i, (N + 1) * i + 1] = 1;
                matr[(N + 1) * i + N, (N + 1) * i + N] = 1;
                d[(N + 1) * i + N] = Math.Pow(h * i, 2) - 1;
                for (int j = 1; j <= N - 1; j++)
                {
                    matr[(N + 1) * i + j, (N + 1) * i - N - 1 + j] = 1;
                    matr[(N + 1) * i + j, (N + 1) * i + N + 1 + j] = 1;
                    matr[(N + 1) * i + j, (N + 1) * i + j] = -4;
                    matr[(N + 1) * i + j, (N + 1) * i - 1 + j] = 1;
                    matr[(N + 1) * i + j, (N + 1) * i + 1 + j] = 1;
                }
            }
            double[] xP = Prost(matr, d);
            double[] xZ = Zeydel(matr, d);
            double[] xPU = ProstUpgrade(matr, d);
            for (int i = 0; i < N + 1; i++)
            {
                for (int j = 0; j < N + 1; j++)
                {
                    U[i, j] = xP[(N + 1) * i + j];
                    dt[i] += Math.Abs(U[i, j] - f(h * i, h * j));
                    U2[i, j] = xZ[(N + 1) * i + j];
                    dt2[i] += Math.Abs(U2[i, j] - f(h * i, h * j));
                    U3[i, j] = xPU[(N + 1) * i + j];
                    dt3[i] += Math.Abs(U3[i, j] - f(h * i, h * j));
                }
            }
            Console.Write(Pcount + "\t" + Zcount + "\t" + PUcount);
            Console.WriteLine();
        }

        public static double[] Prost(double[,] matr, double[] d)
        {
            bool flag = false;
            double[,] B = new double[d.Length, d.Length];
            double[] dup = new double[d.Length];
            double[] x = new double[d.Length];
            double[] x1 = new double[d.Length];
            double e = 0.001;
            double max = 0;
            for (int i = 0; i < d.Length; i++)
            {
                for (int j = 0; j < d.Length; j++) if (i != j) if (matr[i,i] != 0) B[i, j] = -matr[i, j] / matr[i, i];
                if (matr[i, i] != 0) dup[i] = d[i] / matr[i, i];
            }
            for (int i = 0; i < d.Length; i++)
            {
                x[i] = dup[i];
            }
            while (flag == false)
            {
                for (int i = 0; i < d.Length; i++)
                {
                    for (int j = 0; j < d.Length; j++)
                    {
                        x1[i] += B[i, j] * x[j];
                    }
                    x1[i] += dup[i];
                }
                for (int i = 0; i < d.Length; i++)
                {
                    if (Math.Abs(x[i] - x1[i]) > max) max = Math.Abs(x[i] - x1[i]);
                }
                if (max < e) flag = true;
                else
                {
                    for (int i = 0; i < d.Length; i++)
                    {
                        x[i] = x1[i];
                        x1[i] = 0;
                    }
                }
                max = 0;
                Pcount++;
            }
            return x;
        }
        public static double[] Zeydel(double[,] matr, double[] d)
        {
            bool flag = false;
            double[,] B = new double[d.Length, d.Length];
            double[] dup = new double[d.Length];
            double[] x = new double[d.Length];
            double[] x1 = new double[d.Length];
            double[] x2 = new double[d.Length];
            double e = 0.001;
            double max = 0;
            for (int i = 0; i < d.Length; i++)
            {
                for (int j = 0; j < d.Length; j++) if (i != j) if (matr[i, i] != 0) B[i, j] = -matr[i, j] / matr[i, i];
                if (matr[i, i] != 0) dup[i] = d[i] / matr[i, i];
            }
            for (int i = 0; i < d.Length; i++)
            {
                x[i] = dup[i];
                x2[i] = dup[i];
            }
            while (flag == false)
            {
                for (int i = 0; i < d.Length; i++)
                {
                    for (int j = 0; j < d.Length; j++)
                    {
                        x1[i] += B[i, j] * x2[j];
                    }
                    x1[i] += dup[i];
                    x2[i] = x1[i];
                }
                for (int i = 0; i < d.Length; i++)
                {
                    if (Math.Abs(x[i] - x1[i]) > max) max = Math.Abs(x[i] - x1[i]);
                }
                if (max < e) flag = true;
                else
                {
                    for (int i = 0; i < d.Length; i++)
                    {
                        x[i] = x1[i];
                        x2[i] = x1[i];
                        x1[i] = 0;
                    }
                }
                max = 0;
                Zcount++;
            }
            return x;
        }
        public static double[] ProstUpgrade(double[,] matr, double[] d)
        {
            bool flag = false;
            double[,] B = new double[d.Length, d.Length];
            double[] dup = new double[d.Length];
            double[] x = new double[d.Length];
            double[] x1 = new double[d.Length];
            double[] x2 = new double[d.Length];
            double e = 0.001;
            double w = 1.02865;
            double max = 0;
            for (int i = 0; i < d.Length; i++)
            {
                for (int j = 0; j < d.Length; j++) if (i != j) if (matr[i, i] != 0) B[i, j] = -matr[i, j] / matr[i, i];
                if (matr[i, i] != 0) dup[i] = d[i] / matr[i, i];
            }
            for (int i = 0; i < d.Length; i++)
            {
                x[i] = dup[i];
                x2[i] = dup[i];
            }
            while (flag == false)
            {
                for (int i = 0; i < d.Length; i++)
                {
                    for (int j = 0; j < d.Length; j++)
                    {
                        x1[i] += B[i, j] * x2[j];
                    }
                    x1[i] += dup[i];
                }
                for (int i = 0; i < d.Length; i++)
                {
                    x2[i] = x1[i] + (w - 1) * (x1[i] - x[i]);
                }
                for (int i = 0; i < d.Length; i++)
                {
                    if (Math.Abs(x[i] - x2[i]) > max) max = Math.Abs(x[i] - x1[i]);
                }
                if (max < e) flag = true;
                else
                {
                    for (int i = 0; i < d.Length; i++)
                    {
                        x[i] = x2[i];
                        x1[i] = 0;
                    }
                }
                max = 0;
                PUcount++;
            }
            return x;
        }
        private double f(double x, double y) // Точное значение функции
        {
            return Math.Pow(x,2) - Math.Pow(y,2);
        }
        private void chart1_Paint(object sender, PaintEventArgs e)
        {
            chart1.Series[1].Points.Clear();
            chart1.Series[2].Points.Clear();
            chart1.Series[3].Points.Clear();
            if (M >= 0)
            {
                textBox1.Text = (Math.Round(dt[M] / (N + 1), 15)).ToString();
                textBox2.Text = (Math.Round(dt2[M] / (N + 1), 15)).ToString();
                textBox3.Text = (Math.Round(dt3[M] / (N + 1), 15)).ToString();
            }
            if (radioButton3.Checked) paint();
            else if (radioButton2.Checked)
            {
                for (int i = 0; i <= N; i++)
                {
                    chart1.Series[1].Points.AddXY(i * l / N, U[i, M]);
                    chart1.Series[2].Points.AddXY(i * l / N, U2[i, M]);
                    chart1.Series[3].Points.AddXY(i * l / N, U3[i, M]);
                }
                paintY();
            }
            else if (radioButton1.Checked)
            {
                for (int j = 0; j <= N; j++)
                {
                    chart1.Series[1].Points.AddXY(j * l / N, U[M, j]);
                    chart1.Series[2].Points.AddXY(j * l / N, U2[M, j]);
                    chart1.Series[3].Points.AddXY(j * l / N, U3[M, j]);
                }
                paintX();
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {

            U = new double[K + 1, N + 1];
            dt = new double[K + 1];
            CRS();
        }

        private void trackBar1_Scroll(object sender, EventArgs e)
        {
            M = trackBar1.Value;
            textBox4.Text = (M / Convert.ToDouble(K)).ToString();
        }

        private void radioButton3_CheckedChanged(object sender, EventArgs e)
        {
            trackBar1.Minimum = -K;
            trackBar1.Value = 0;
            M = 0;
            textBox4.Text = (M / Convert.ToDouble(K)).ToString();
            label2.Text = "Значение U(x,y):";
        }

        private void radioButton2_CheckedChanged(object sender, EventArgs e)
        {
            trackBar1.Minimum = 0;
            trackBar1.Value = 0;
            M = 0;
            textBox4.Text = (0).ToString();
            label2.Text = "Значение Y:";
        }

        private void radioButton1_CheckedChanged(object sender, EventArgs e)
        {
            trackBar1.Minimum = 0;
            trackBar1.Value = 0;
            M = 0;
            textBox4.Text = (0).ToString();
            label2.Text = "Значение X:";
        }
    }
}
