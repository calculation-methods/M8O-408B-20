using System;
using System.Windows.Forms;
using System.Collections.Generic;

namespace NML5
{
    public partial class Form1 : Form
    {
        
        int M;
        int N = 100;
        int K = 100;
        double l = Math.PI; 
        double[,] U = new double[101,101];
        double[,] U2 = new double[101, 101];
        double[,] U3 = new double[101, 101];
        double[] dt = new double[101];
        double[] dt2 = new double[101];
        double[] dt3 = new double[101];
        double a = 1;
        double sig = 0.9;
        public Form1()
        {
            InitializeComponent();
        }
        private void paint()
        {
            double h = l / N;
            double tau = Math.Sqrt(sig * Math.Pow(h, 2) / Math.Pow(a, 2));
            double u, x;
            chart1.Series[0].Points.Clear();
            chart1.ChartAreas[0].AxisY.LabelStyle.Format = "f2";
            chart1.ChartAreas[0].AxisX.LabelStyle.Format = "f2";
            chart1.ChartAreas[0].AxisX.Maximum = l;
            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisY.Maximum = 2; 
            chart1.ChartAreas[0].AxisY.Minimum = -2;
            for (double i = 0; i <= Math.PI; i += h)
            {
                x = i;
                u = Math.Sin(x - a * M * tau) + Math.Cos(x + a * M * tau); 
                chart1.Series[0].Points.AddXY(x, u);
            }
        }
        private void Yav(bool flag, int apr) //Явная схема
        {
            double h = l / N;
            double tau = Math.Sqrt(sig * Math.Pow(h, 2) / Math.Pow(a,2));
            double[,] u = new double[K + 1, N + 1];
            for (int j = 0; j <= N; j++)
            {
                u[0, j] = Math.Sin(j * h) + Math.Cos(j * h);
                if (flag)
                    u[1, j] = Math.Sin(j * h) + Math.Cos(j * h) - a * (Math.Sin(j * h) + Math.Cos(j * h)) * tau;
                else
                    u[1, j] = Math.Sin(j * h) + Math.Cos(j * h) - a * (Math.Sin(j * h) + Math.Cos(j * h)) * tau + Math.Pow(a, 2) * (-Math.Sin(j * h) - Math.Cos(j * h)) * Math.Pow(tau, 2) / 2;
            }
            for (int k = 1; k <= K - 1; k ++)
            {
                for (int j = 1; j <= N-1; j++)
                {
                    u[k + 1, j] = u[k, j + 1] * sig + u[k, j] * (-2 * sig + 2) + u[k, j - 1] * sig - u[k - 1, j];
                }
                if (apr == 0)
                {
                    u[k + 1, 0] = u[k + 1, 1] / (h + 1);
                    u[k + 1, N] = u[k + 1, N - 1] / (1 - h);
                }
                else if (apr == 1) 
                {
                    u[k + 1, 0] = (-4 * u[k + 1, 1] + u[k + 1, 2]) / (-2 * h - 3);
                    u[k + 1, N] = (4 * u[k + 1, N - 1] - u[k + 1, N - 2]) / (-2 * h + 3);
                }
                else if (apr == 2) 
                {
                    u[k + 1, 0] = (u[k + 1, 1] + 1 / (2 * sig) * (2 * u[k, 0] - u[k - 1, 0])) / (h + 1 + 1 / (2 * sig));
                    u[k + 1, N] = (u[k + 1, N - 1] + 1 / (2 * sig) * (2 * u[k, N] - u[k - 1, N])) / (-h + 1 + 1 / (2 * sig));
                }
            }
            if (apr == 0)
            {
                for (int k = 0; k <= K; k++)
                {
                    for (int j = 0; j <= N; j++)
                    {
                        U[k, j] = u[k, j];
                        dt[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                    }
                }
            }
            else if (apr == 1)
            {
                for (int k = 0; k <= K; k++)
                {
                    for (int j = 0; j <= N; j++)
                    {
                        U2[k, j] = u[k, j];
                        dt2[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                    }
                }
            }
            else if (apr == 2)
            {
                for (int k = 0; k <= K; k++)
                {
                    for (int j = 0; j <= N; j++)
                    {
                        U3[k, j] = u[k, j];
                        dt3[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                    }
                }
            }
        }
        
        private void NeYav(bool flag, int apr) //Неявная схема 2Т1П
        {
            double h = l / N;
            double tau = Math.Sqrt(sig * Math.Pow(h, 2) / Math.Pow(this.a, 2));
            double[,] u = new double[K + 1, N + 1];
            double[] b = new double[N + 1];
            double[] a = new double[N];
            double[] c = new double[N];
            double[] d = new double[N + 1];
            double[] x;
            double p;
            for (int j = 0; j <= N; j++)
            {
                u[0, j] = Math.Sin(j * h) + Math.Cos(j * h);
                if (flag)
                    u[1, j] = Math.Sin(j * h) + Math.Cos(j * h) - this.a * (Math.Sin(j * h) + Math.Cos(j * h)) * tau;
                else
                    u[1, j] = Math.Sin(j * h) + Math.Cos(j * h) - this.a * (Math.Sin(j * h) + Math.Cos(j * h)) * tau + Math.Pow(this.a, 2) * (-Math.Sin(j * h) - Math.Cos(j * h)) * Math.Pow(tau, 2) / 2;
            }
            for (int k = 1; k <= K - 1; k++)
            {
                for (int j = 0; j <= N - 1; j++)
                {
                    a[j] = -sig; b[j] = 1 + 2 * sig; c[j] = -sig; d[j] = 2 * u[k, j] - u[k - 1, j];
                }
                if (apr == 0)
                {
                    b[0] = -1 / h - 1; c[0] = 1/h; d[0] = 0;
                    a[N - 1] = -1 / h; b[N] = 1 / h - 1; d[N] = 0; 
                }
                else if (apr == 1)
                {
                    p = 1 / (2 * h * sig);
                    b[0] = -3 / (2 * h) -1 - p * a[0]; c[0] = 4 / (2 * h) - p * b[1]; d[0] = - p * d[1]; 
                    a[N - 1] = -4 / (2 * h) + p * b[N - 1]; b[N] = 3 / (2 * h) - 1 + p * c[N - 1]; d[N] =  p * d[N - 1]; 
                }
                else if (apr == 2)
                {
                    b[0] = 1 + 1 / (2 * sig) + h; c[0] = -1; d[0] = 1 / (2 * sig) * (2 * u[k, 0] - u[k - 1, 0]);
                    a[N - 1] = -1; b[N] = 1 + 1 / (2 * sig) - h; d[N] = 1 / (2 * sig) * (2 * u[k, N] - u[k - 1, N]); 
                }
                x = Progon(a,b,c,d).ToArray();
                for (int j = 0; j <= N; j++)
                {
                    u[k + 1, j] = x[j];
                }
            }
            if (apr == 0)
            {
                for (int k = 0; k <= K; k++)
                {
                    for (int j = 0; j <= N; j++)
                    {
                        U[k, j] = u[k, j];
                        dt[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                    }
                }
            }
            else if (apr == 1)
            {
                for (int k = 0; k <= K; k++)
                {
                    for (int j = 0; j <= N; j++)
                    {
                        U2[k, j] = u[k, j];
                        dt2[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                    }
                }
            }
            else if (apr == 2)
            {
                for (int k = 0; k <= K; k++)
                {
                    for (int j = 0; j <= N; j++)
                    {
                        U3[k, j] = u[k, j];
                        dt3[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                    }
                }
            }
        }
        static List<Double> Progon(double[] arr, double[] b, double[] crr, double[] d) // Метод прогонки
        {
            double[] a = new double[arr.Length + 1];
            double[] c = new double[arr.Length + 1];
            a[0] = 0; c[crr.Length] = 0;
            for (int i = 1; i < a.Length; i++) a[i] = arr[i - 1];
            for (int i = 0; i < a.Length - 1; i++) c[i] = crr[i];
            List<Double> roots = new List<double>();
            List<Double> P = new List<double>();
            List<Double> Q = new List<double>();
            P.Add(-c[0] / b[0]);
            Q.Add(d[0] / b[0]);
            for (int i = 1; i < a.Length; i++)
            {
                P.Add(-c[i] / (b[i] + a[i] * P[(i - 1)]));
                Q.Add((d[i] - a[i] * Q[(i - 1)]) / (b[i] + a[i] * P[(i - 1)]));
            }
            P.Reverse();
            Q.Reverse();
            roots.Add(Q[0]);
            for (int i = 1; i < a.Length; i++)
            {
                roots.Add(P[i] * roots[i - 1] + Q[i]);
            }
            roots.Reverse();
            return roots;
        }
        private double f(double x, double t) // Точное значение функции
        {
            return Math.Sin(x - a * t) + Math.Cos(x + a * t);
        }
        private void chart1_Paint(object sender, PaintEventArgs e)
        {
            paint();
            chart1.Series[1].Points.Clear();
            textBox1.Text = (Math.Round(dt[M] / (N + 1),15)).ToString();
            if (checkBox1.Checked == true)
            {
                for (int j = 0; j <= N; j++)
                {
                    chart1.Series[1].Points.AddXY(j * l / N, U[M, j]);
                }
            }
            chart1.Series[2].Points.Clear();
            textBox2.Text = (Math.Round(dt2[M] / (N + 1), 15)).ToString();
            if (checkBox2.Checked == true)
            {
                for (int j = 0; j <= N; j++)
                {
                    chart1.Series[2].Points.AddXY(j * l / N, U2[M, j]);
                }
            }
            chart1.Series[3].Points.Clear();
            textBox3.Text = (Math.Round(dt3[M] / (N + 1), 15)).ToString();
            if (checkBox3.Checked == true)
            {
                for (int j = 0; j <= N; j++)
                {
                    chart1.Series[3].Points.AddXY(j * l / N, U3[M, j]);
                }
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {
            bool flag = true;
            if (radioButton1.Checked)
            {
                if (radioButton3.Checked) flag = true;
                else if (radioButton4.Checked) flag = false;
                U = new double[K + 1, N + 1];
                dt = new double[K + 1];
                U3 = new double[K + 1, N + 1];
                dt3 = new double[K + 1];
                U2 = new double[K + 1, N + 1];
                dt2 = new double[K + 1];
                Yav(flag, 0);
                Yav(flag, 1);
                Yav(flag, 2);
            }
            else if (radioButton2.Checked)
            {
                if (radioButton3.Checked) flag = true;
                else if (radioButton4.Checked) flag = false;
                U = new double[K + 1, N + 1];
                dt = new double[K + 1];
                U3 = new double[K + 1, N + 1];
                dt3 = new double[K + 1];
                U2 = new double[K + 1, N + 1];
                dt2 = new double[K + 1];
                NeYav(flag, 0);
                NeYav(flag, 1);
                NeYav(flag, 2);
            }
        }

        private void trackBar1_Scroll(object sender, EventArgs e)
        {
            M = trackBar1.Value;
            textBox4.Text = M.ToString();
        }
    }
}
