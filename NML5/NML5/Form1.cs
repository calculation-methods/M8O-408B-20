using System;
using System.Windows.Forms;
using System.Collections.Generic;

namespace NML5
{
    public partial class Form1 : Form
    {
        
        int M;
        int N = 10;
        int K = 100;
        double l = 1; //2-й вариант
        //double l = Math.PI; //4-й вариант
        double[,] U = new double[101,101];
        double[,] U2 = new double[101, 101];
        double[,] U3 = new double[101, 101];
        double[] dt = new double[101];
        double[] dt2 = new double[101];
        double[] dt3 = new double[101];
        double a = 1;
        double sig = 0.5;
        public Form1()
        {
            InitializeComponent();
        }
        private void paint()
        {
            double h = l / N;
            double tau = sig * Math.Pow(h, 2) / a;
            double u, x;
            chart1.Series[0].Points.Clear();
            chart1.ChartAreas[0].AxisY.LabelStyle.Format = "f2";
            chart1.ChartAreas[0].AxisX.LabelStyle.Format = "f2";
            chart1.ChartAreas[0].AxisX.Maximum = l;
            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisY.Maximum = 1.7; //2-й вариант
            //chart1.ChartAreas[0].AxisY.Maximum = 1; //4-й вариант
            chart1.ChartAreas[0].AxisY.Minimum = 0;
            for (double i = 0; i <= Math.PI; i += h)
            {
                x = i;
                //u = Math.Exp(-a * M * tau) * Math.Sin(x); //4-й вариант
                u = x + Math.Exp(-Math.Pow(Math.PI, 2) * a * M * tau) * Math.Sin(Math.PI * x); //2-й вариант
                chart1.Series[0].Points.AddXY(x, u);
            }
        }
        private void Yav2D() //Явная схема 2Т1П
        {
            double h = l / N;
            double tau = sig * Math.Pow(h, 2) / a;
            double[,] u = new double[K + 1, N + 1];
            for (int j = 0; j <= N; j++) u[0, j] = j * h + Math.Sin(Math.PI * j * h); //2-й вариант
            //for (int j = 0; j <= N; j++) u[0, j] = Math.Sin(j * h); //4-й вариант
            for (int k = 0; k <= K - 1; k ++)
            {
                for (int j = 1; j <= N-1; j++)
                {
                    u[k + 1, j] = sig * u[k, j + 1] + (1 - 2 * sig) * u[k, j] + sig * u[k, j - 1];
                }
                u[k + 1, 0] = 0; //2-й вариант
                u[k + 1, N] = 1; //2-й вариант
                //u[k + 1, 0] = -(h * Math.Exp(-a * (k + 1) * tau) - u[k + 1, 1]); //4-й вариант
                //u[k + 1, N] = -h * Math.Exp(-a * (k + 1) * tau) + u[k + 1, N - 1]; //4-й вариант
            }
            for (int k = 0; k <= K; k++)
            {
                for (int j = 0; j <= N; j++)
                {
                    U[k, j] = u[k, j];
                    dt[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                }
            }
        }
        private void Yav3D() //Явная схема 3Т3П
        {
            double h = l / N;
            double tau = sig * Math.Pow(h, 2) / a;
            double[,] u = new double[K + 1, N + 1];
            for (int j = 0; j <= N; j++) u[0, j] = j * h + Math.Sin(Math.PI * j * h); //2-й вариант
            //for (int j = 0; j <= N; j++) u[0, j] = Math.Sin(j * h); //4-й вариант
            for (int k = 0; k <= K - 1; k++)
            {
                for (int j = 1; j <= N - 1; j++)
                {
                    u[k + 1, j] = sig * u[k, j + 1] + (1 - 2 * sig) * u[k, j] + sig * u[k, j - 1];
                }
                u[k + 1, 0] = 0; //2-й вариант
                u[k + 1, N] = 1; //2-й вариант
                //u[k + 1, 0] = (2 * h * Math.Exp(-a * (k + 1) * tau) - 4 * u[k + 1, 1] + u[k + 1, 2]) / (-3); //4-й вариант
                //u[k + 1, N] = (2 * h * (-Math.Exp(-a * (k + 1) * tau)) + 4 * u[k + 1, N - 1] - u[k + 1, N - 2]) / 3; //4-й вариант
            }
            for (int k = 0; k <= K; k++)
            {
                for (int j = 0; j <= N; j++)
                {
                    U2[k, j] = u[k, j];
                    dt2[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                }
            }
        }
        private void Yav2B2T() //Явная схема 2Т2П
        {
            double h = l / N;
            double tau = sig * Math.Pow(h, 2) / a;
            double[,] u = new double[K + 1, N + 1];
            for (int j = 0; j <= N; j++) u[0, j] = j * h + Math.Sin(Math.PI * j * h); //2-й вариант
            //for (int j = 0; j <= N; j++) u[0, j] = Math.Sin(j * h); //4-й вариант
            for (int k = 0; k <= K - 1; k++)
            {
                for (int j = 1; j <= N - 1; j++)
                {
                    u[k + 1, j] = sig * u[k, j + 1] + (1 - 2 * sig) * u[k, j] + sig * u[k, j - 1];
                }
                u[k + 1, 0] = 0; //2-й вариант
                u[k + 1, N] = 1; //2-й вариант
                //u[k + 1, 0] = (Math.Exp(-a * (k + 1) * tau) - u[k + 1, 1] / h - h * u[k, 0] / (2 * a * tau)) / (-1 / h - h / (2 * a * tau)); //4-й вариант
                //u[k + 1, N] = (-Math.Exp(-a * (k + 1) * tau) + u[k + 1, N - 1] / h + h * u[k, N] / (2 * a * tau)) / (1 / h + h / (2 * a * tau)); //4-й вариант
            }
            for (int k = 0; k <= K; k++)
            {
                for (int j = 0; j <= N; j++)
                {
                    U3[k, j] = u[k, j];
                    dt3[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                }
            }
        }
        private void NeYav2D() //Неявная схема 2Т1П
        {
            double h = l / N;
            double tau = sig * Math.Pow(h, 2) / this.a;
            double[,] u = new double[K + 1, N + 1];
            double[] b = new double[N + 1];
            double[] a = new double[N];
            double[] c = new double[N];
            double[] d = new double[N + 1];
            double[] x;
            for (int j = 0; j <= N; j++) u[0, j] = j * h + Math.Sin(Math.PI * j * h); //2-й вариант
            //for (int j = 0; j <= N; j++) u[0, j] = Math.Sin(j * h); //4-й вариант
            for (int k = 0; k <= K - 1; k++)
            {
                for (int j = 0; j <= N - 1; j++)
                {
                    a[j] = sig; b[j] = -(1 + 2 * sig); c[j] = sig; d[j] = -u[k, j];
                }
                b[0] = 1; c[0] = 0; d[0] = 0;  //2-й вариант
                a[N - 1] = 0; b[N] = 1; d[N] = 1; //2-й вариант
                //b[0] = -1/h; c[0] = 1/h; d[0] = Math.Exp(-this.a * (k + 1) * tau); //4-й вариант
                //a[N - 1] = -1/h; b[N] = 1/h; d[N] = (-Math.Exp(-this.a * (k + 1) * tau)); //4-й вариант
                x = Progon(a,b,c,d).ToArray();
                for (int j = 0; j <= N; j++)
                {
                    u[k + 1, j] = x[j];
                }
            }
            for (int k = 0; k <= K; k++)
            {
                for (int j = 0; j <= N; j++)
                {
                    U[k, j] = u[k, j];
                    dt[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                }
            }
        }
        private void NeYav3D() //Неявная схема 3Т2П
        {
            double h = l / N;
            double tau = sig * Math.Pow(h, 2) / this.a;
            double[,] u = new double[K + 1, N + 1];
            double[] b = new double[N + 1];
            double[] a = new double[N];
            double[] c = new double[N];
            double[] d = new double[N + 1];
            double[] x;
            double p, p1;
            for (int j = 0; j <= N; j++) u[0, j] = j * h + Math.Sin(Math.PI * j * h);
            //for (int j = 0; j <= N; j++) u[0, j] = Math.Sin(j * h); //4-й вариант
            for (int k = 0; k <= K - 1; k++)
            {
                for (int j = 0; j <= N - 1; j++)
                {
                    a[j] = sig; b[j] = -(1 + 2 * sig); c[j] = sig; d[j] = -u[k, j];
                }
                b[0] = 1; c[0] = 0; d[0] = 0; //2-й вариант
                a[N - 1] = 0; b[N] = 1; d[N] = 1; //2-й вариант
                //p = (-1 / (2 * h)) / sig; p1 = (1 / (2 * h)) / sig;
                //b[0] = -3 / (2 * h) - p * a[0]; c[0] = 4 / (2 * h) - p * b[1]; d[0] = Math.Exp(-this.a * (k + 1) * tau) - p * d[1]; //4-й вариант
                //a[N - 1] = -4 / (2 * h) - p1 * b[N - 1]; b[N] = 3 / (2 * h) - p1 * c[N - 1]; d[N] = (-Math.Exp(-this.a * (k + 1) * tau)) - p1 * d[N - 1]; //4-й вариант
                x = Progon(a, b, c, d).ToArray(); 
                for (int j = 0; j <= N; j++)
                {
                    u[k + 1, j] = x[j];
                }
            }
            for (int k = 0; k <= K; k++)
            {
                for (int j = 0; j <= N; j++)
                {
                    U2[k, j] = u[k, j];
                    dt2[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                }
            }
        }
        private void NeYav2B2T() //Неявная схема 2Т2П
        {
            double h = l / N;
            double tau = sig * Math.Pow(h, 2) / this.a;
            double[,] u = new double[K + 1, N + 1];
            double[] b = new double[N + 1];
            double[] a = new double[N];
            double[] c = new double[N];
            double[] d = new double[N + 1];
            double[] x;
            for (int j = 0; j <= N; j++) u[0, j] = j * h + Math.Sin(Math.PI * j * h); //2-й вариант
            //for (int j = 0; j <= N; j++) u[0, j] = Math.Sin(j * h); //4-й вариант
            for (int k = 0; k <= K - 1; k++)
            {
                for (int j = 0; j <= N - 1; j++)
                {
                    a[j] = -this.a / Math.Pow(h, 2); b[j] = 2 * this.a / Math.Pow(h, 2) + 1 / tau; c[j] = -this.a / Math.Pow(h, 2); d[j] = 1 / tau * u[k, j];
                }
                b[0] = 1; c[0] = 0; d[0] = 0; //2-й вариант
                a[N - 1] = 0; b[N] = 1; d[N] = 1; //2-й вариант
                //b[0] = 2 * this.a / h + h / tau; c[0] = -2 * this.a / h; d[0] = h / tau * u[k, 0] - Math.Exp(-this.a * (k + 1) * tau) * 2 * this.a; //4-й вариант
                //a[N - 1] = -2 * this.a / h; b[N] = 2 * this.a / h + h / tau; d[N] = h / tau * u[k, N] + (-Math.Exp(-this.a * (k + 1) * tau)) * 2 * this.a; //4-й вариант
                x = Progon(a, b, c, d).ToArray();
                for (int j = 0; j <= N; j++)
                {
                    u[k + 1, j] = x[j];
                }
            }
            for (int k = 0; k <= K; k++)
            {
                for (int j = 0; j <= N; j++)
                {
                    U3[k, j] = u[k, j];
                    dt3[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                }
            }
        }
        private void KN2D() //Схема Кранка-Николсона 2Т1П
        {
            double h = l / N;
            double tau = sig * Math.Pow(h, 2) / this.a;
            double[,] u = new double[K + 1, N + 1];
            double[] b = new double[N + 1];
            double[] a = new double[N];
            double[] c = new double[N];
            double[] d = new double[N + 1];
            double[] x;
            double r = this.a * tau / (h * h);
            for (int j = 0; j <= N; j++) u[0, j] = j * h + Math.Sin(Math.PI * j * h); //2-й вариант
            //for (int j = 0; j <= N; j++) u[0, j] = Math.Sin(j * h); //4-й вариант
            for (int k = 0; k <= K - 1; k++)
            {
                for (int j = 1; j <= N - 1; j++)
                {
                    a[j] = -r/2; b[j] = r+1; c[j] = -r/2; d[j] = r / 2 * (u[k,j - 1] + u[k,j + 1]) + u[k,j] * (1 - r);
                }
                b[0] = 1; c[0] = 0; d[0] = 0; a[0] = -r / 2; //2-й вариант
                a[N - 1] = 0; b[N] = 1; d[N] = 1; c[N - 1] = -r / 2;//2-й вариант
                //b[0] = -1 / h; c[0] = 1 / h; d[0] = Math.Exp(-this.a * (k + 1) * tau); a[0] = -r / 2; //4-й вариант
                //a[N - 1] = -1 / h; b[N] = 1 / h; d[N] = (-Math.Exp(-this.a * (k + 1) * tau)); c[N - 1] = -r / 2; //4-й вариант
                x = Progon(a, b, c, d).ToArray();
                for (int j = 0; j <= N; j++)
                {
                    u[k + 1, j] = x[j];
                }
            }
            for (int k = 0; k <= K; k++)
            {
                for (int j = 0; j <= N; j++)
                {
                    U[k, j] = u[k, j];
                    dt[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                }
            }
        }
        private void KN3D() //Схема Кранка-Николсона 3Т2П
        {
            double h = l / N;
            double tau = sig * Math.Pow(h, 2) / this.a;
            double[,] u = new double[K + 1, N + 1];
            double[] b = new double[N + 1];
            double[] a = new double[N];
            double[] c = new double[N];
            double[] d = new double[N + 1];
            double[] x;
            double koef = 1 / (2 * h * sig);
            double r = this.a * tau / (h * h);
            for (int j = 0; j <= N; j++) u[0, j] = j * h + Math.Sin(Math.PI * j * h); //2-й вариант
            //for (int j = 0; j <= N; j++) u[0, j] = Math.Sin(j * h); //4-й вариант
            for (int k = 0; k <= K - 1; k++)
            {
                for (int j = 1; j <= N - 1; j++)
                {
                    a[j] = -r / 2; b[j] = r + 1; c[j] = -r / 2; d[j] = r / 2 * (u[k, j - 1] + u[k, j + 1]) + u[k, j] * (1 - r);
                }
                b[0] = 1; c[0] = 0; d[0] = 0; a[0] = -r / 2; //2-й вариант
                a[N - 1] = 0; b[N] = 1; d[N] = 1; c[N - 1] = -r / 2;//2-й вариант
                //b[0] = -1 / h; c[0] = 1 / h - koef; d[0] = Math.Exp(-this.a * (k + 1) * tau) - koef * u[k, 1]; a[0] = -r / 2; //4-й вариант
                //a[N - 1] = -1 / h + koef; b[N] = 1 / h; d[N] = -Math.Exp(-this.a * (k + 1) * tau) + koef * u[k,N - 1]; c[N - 1] = -r / 2; //4-й вариант
                x = Progon(a, b, c, d).ToArray();
                for (int j = 0; j <= N; j++)
                {
                    u[k + 1, j] = x[j];
                }
            }
            for (int k = 0; k <= K; k++)
            {
                for (int j = 0; j <= N; j++)
                {
                    U2[k, j] = u[k, j];
                    dt2[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
                }
            }
        }
        private void KN2B2T() //Схема Кранка-Николсона 2Т2П
        {
            double h = l / N;
            double tau = sig * Math.Pow(h, 2) / this.a;
            double[,] u = new double[K + 1, N + 1];
            double[] b = new double[N + 1];
            double[] a = new double[N];
            double[] c = new double[N];
            double[] d = new double[N + 1];
            double[] x;
            double r = this.a * tau / (h * h);
            for (int j = 0; j <= N; j++) u[0, j] = j * h + Math.Sin(Math.PI * j * h); //2-й вариант
            //for (int j = 0; j <= N; j++) u[0, j] = Math.Sin(j * h); //4-й вариант
            for (int k = 0; k <= K - 1; k++)
            {
                for (int j = 1; j <= N - 1; j++)
                {
                    a[j] = -r / 2; b[j] = r + 1; c[j] = -r / 2; d[j] = r / 2 * (u[k, j - 1] + u[k, j + 1]) + u[k, j] * (1 - r);
                }
                b[0] = 1; c[0] = 0; d[0] = 0; a[0] = -r / 2; //2-й вариант
                a[N - 1] = 0; b[N] = 1; d[N] = 1; c[N - 1] = -r / 2; //2-й вариант
                //b[0] = 2 * this.a / h + h / tau; c[0] = -2 * this.a / h; d[0] = h / tau * u[k,0] - 2 * this.a * Math.Exp(-this.a * (k + 1) * tau); a[0] = -r / 2; //4-й вариант
                //a[N - 1] = -2 * this.a / h; b[N] = 2 * this.a / h + h / tau; d[N] = h / tau * u[k,N] - 2 * this.a * Math.Exp(-this.a * (k + 1) * tau); c[N - 1] = -r / 2; //4-й вариант
                x = Progon(a, b, c, d).ToArray();
                for (int j = 0; j <= N; j++)
                {
                    u[k + 1, j] = x[j];
                }
            }
            for (int k = 0; k <= K; k++)
            {
                for (int j = 0; j <= N; j++)
                {
                    U3[k, j] = u[k, j];
                    dt3[k] += Math.Abs(u[k, j] - f(h * j, tau * k));
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
            return x + Math.Exp(-Math.Pow(Math.PI, 2) * a * t) * Math.Sin(Math.PI * x); //2-й вариант
            //return Math.Exp(-a * t) * Math.Sin(x); //4-й вариант
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
            if (radioButton1.Checked)
            {
                U = new double[K + 1, N + 1];
                dt = new double[K + 1];
                U3 = new double[K + 1, N + 1];
                dt3 = new double[K + 1];
                U2 = new double[K + 1, N + 1];
                dt2 = new double[K + 1];
                Yav2D();
                Yav3D();
                Yav2B2T();
            }
            else if (radioButton2.Checked)
            {
                U = new double[K + 1, N + 1];
                dt = new double[K + 1];
                U3 = new double[K + 1, N + 1];
                dt3 = new double[K + 1];
                U2 = new double[K + 1, N + 1];
                dt2 = new double[K + 1];
                NeYav2D();
                NeYav3D();
                NeYav2B2T();
            }
            else if (radioButton3.Checked)
            {
                U = new double[K + 1, N + 1];
                dt = new double[K + 1];
                U3 = new double[K + 1, N + 1];
                dt3 = new double[K + 1];
                U2 = new double[K + 1, N + 1];
                dt2 = new double[K + 1];
                KN2D();
                KN3D();
                KN2B2T();
            }
        }

        private void trackBar1_Scroll(object sender, EventArgs e)
        {
            M = trackBar1.Value;
            textBox4.Text = M.ToString();
        }
    }
}
