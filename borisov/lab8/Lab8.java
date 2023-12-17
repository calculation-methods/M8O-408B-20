import java.util.ArrayList;
import java.util.Collections;

public class Lab8 {
    public static void main(String[] args) {
        double a = 1;
        int I = 10;
        double hx = Math.PI / (4 * I);
        int J = 20;
        double hy = Math.log(2) / J;
        int T = 20;
        double ht = 1. / T;

        double[][][] u = MPN(a, I, J, T, hx, hy, ht);
        int s = 10;
        int r = 12;

        double max = 0;
        double temp;
        for(int i = 0; i < I + 1; i++){
            temp = Math.abs(u[i][r][s] - U(i * hx, r * hy, s * ht, a));
            if(temp > max){
                max = temp;
            }
        }
        System.out.println("Метод переменных направлений");
        System.out.println(max);

        u = MDS(a, I, J, T, hx, hy, ht);
        max = 0;
        for(int i = 0; i < I + 1; i++){
            temp = Math.abs(u[i][r][s] - U(i * hx, r * hy, s * ht, a));
            if(temp > max){
                max = temp;
            }
        }
        System.out.println("\nМетод дробных шагов");
        System.out.println(max);
    }

    private static double[][][] MPN(double a, int I, int J, int T, double hx, double hy, double ht){
        double [][][] u = new double[I + 1][J + 1][T + 1];
        double [][] u_ = new double[I + 1][J + 1];
        double[] a_arr;
        double[] b_arr;
        double[] c_arr;
        double[] d_arr;

        for (int i = 0; i < I + 1; i++) {
            for (int j = 0; j < J + 1; j++) {
                u[i][j][0] = Math.cos(2 * i * hx) * Math.cosh(j * hy);
            }
        }

        for (int k = 0; k < T + 1; k++) {
            for (int j = 0; j < J + 1; j++) {
                u[0][j][k] = Math.cosh(j * hy) * Math.exp(-3 * a * k * ht);
            }
        }

        for (int k = 0; k < T + 1; k++) {
            for (int j = 0; j < J + 1; j++) {
                u[I][j][k] = 0;
            }
        }

        for (int k = 0; k < T + 1; k++) {
            for (int i = 0; i < I + 1; i++) {
                u[i][0][k] = Math.cos(2 * i * hx) * Math.exp(-3 * a * k * ht);
            }
        }

        double sx = a * ht / (2 * hx * hx);
        double sy = a * ht / (2 * hy * hy);
        int m = 0;

        for(int k = 0; k < T; k++){
            a_arr = new double[I + 1];
            b_arr = new double[I + 1];
            c_arr = new double[I + 1];
            d_arr = new double[I + 1];
            u_ = new double[I + 1][J + 1];
            for(int j = 1; j < J; j++){
                m = 0;
                for(int i = 0; i < I; i++, m++){
                    if(i == 0){
                        a_arr[m] = 0;
                        b_arr[m] = 1;
                        c_arr[m] = 0;
                        d_arr[m] = Math.cosh(j * hy) * Math.exp(-3 * a * (k + 0.5) * ht);
                    }
                    else{
                        a_arr[m] = -sx;
                        b_arr[m] = 1 + 2 * sx;
                        c_arr[m] = -sx;
                        d_arr[m] = u[i][j][k] + sy * (u[i][j + 1][k] - 2 * u[i][j][k] + u[i][j - 1][k]);
                    }
                }

                a_arr[m] = 0;
                b_arr[m] = 1;
                c_arr[m] = 0;
                d_arr[m] = 0;

                ArrayList<Double> result = Progonka(a_arr, b_arr, c_arr, d_arr);
                for(int n = 0; n < u_.length; n++){
                    u_[n][j] = result.get(n);
                }

                if(j == J - 1) {
                    for (int i = 0; i < u_.length; i++) {
                        u_[i][0] = Math.cos(2 * i * hx) * Math.exp(-3 * a * (k + 0.5) * ht);
                    }

                    for (int i = 0; i < u_.length; i++) {
                        u_[i][J] = u_[i][j] + 0.75 * Math.cos(2 * i * hx) * Math.exp(-3 * a * (k + 0.5) * ht);
                    }
                }
            }

            a_arr = new double[J + 1];
            b_arr = new double[J + 1];
            c_arr = new double[J + 1];
            d_arr = new double[J + 1];

            for(int i = 1; i < I; i++){
                m = 0;
                for(int j = 0; j < J; j++, m++){
                    if(j == 0){
                        a_arr[m] = 0;
                        b_arr[m] = 1;
                        c_arr[m] = 0;
                        d_arr[m] = Math.cos(2 * i * hx) * Math.exp(-3 * a * (k + 1) * ht);
                    }
                    else{
                        a_arr[m] = -sy;
                        b_arr[m] = 1 + 2 * sy;
                        c_arr[m] = -sy;
                        d_arr[m] = u_[i][j] + sx * (u_[i + 1][j] - 2 * u_[i][j] + u_[i - 1][j]);
                    }
                }

                a_arr[m] = -1 / hy;
                b_arr[m] = 1 / hy;
                c_arr[m] = 0;
                d_arr[m] = 0.75 * Math.cos(2 * i * hx) * Math.exp(-3 * a * (k + 1) * ht);

                ArrayList<Double> result = Progonka(a_arr, b_arr, c_arr, d_arr);
                for(int j = 0; j < J + 1; j++){
                    u[i][j][k + 1] = result.get(j);
                }
            }
        }
        return u;
    }

    private static double[][][] MDS(double a, int I, int J, int T, double hx, double hy, double ht){
        double [][][] u = new double[I + 1][J + 1][T + 1];
        double [][] u_ = new double[I + 1][J + 1];
        double[] a_arr;
        double[] b_arr;
        double[] c_arr;
        double[] d_arr;

        for (int i = 0; i < I + 1; i++) {
            for (int j = 0; j < J + 1; j++) {
                u[i][j][0] = Math.cos(2 * i * hx) * Math.cosh(j * hy);
            }
        }

        for (int k = 0; k < T + 1; k++) {
            for (int j = 0; j < J + 1; j++) {
                u[0][j][k] = Math.cosh(j * hy) * Math.exp(-3 * a * k * ht);
            }
        }

        for (int k = 0; k < T + 1; k++) {
            for (int j = 0; j < J + 1; j++) {
                u[I][j][k] = 0;
            }
        }

        for (int k = 0; k < T + 1; k++) {
            for (int i = 0; i < I + 1; i++) {
                u[i][0][k] = Math.cos(2 * i * hx) * Math.exp(-3 * a * k * ht);
            }
        }

        double sx = a * ht / (hx * hx);
        double sy = a * ht / (hy * hy);
        int m = 0;

        for(int k = 0; k < T; k++){
            a_arr = new double[I + 1];
            b_arr = new double[I + 1];
            c_arr = new double[I + 1];
            d_arr = new double[I + 1];
            u_ = new double[I + 1][J + 1];
            for(int j = 1; j < J; j++){
                m = 0;
                for(int i = 0; i < I; i++, m++){
                    if(i == 0){
                        a_arr[m] = 0;
                        b_arr[m] = 1;
                        c_arr[m] = 0;
                        d_arr[m] = Math.cosh(j * hy) * Math.exp(-3 * a * (k + 0.5) * ht);
                    }
                    else{
                        a_arr[m] = -sx;
                        b_arr[m] = 1 + 2 * sx;
                        c_arr[m] = -sx;
                        d_arr[m] = u[i][j][k];
                    }
                }

                a_arr[m] = 0;
                b_arr[m] = 1;
                c_arr[m] = 0;
                d_arr[m] = 0;

                ArrayList<Double> result = Progonka(a_arr, b_arr, c_arr, d_arr);
                for(int n = 0; n < u_.length; n++){
                    u_[n][j] = result.get(n);
                }

                if(j == J - 1) {
                    for (int i = 0; i < u_.length; i++) {
                        u_[i][0] = Math.cos(2 * i * hx) * Math.exp(-3 * a * (k + 0.5) * ht);
                    }

                    for (int i = 0; i < u_.length; i++) {
                        u_[i][J] = u_[i][j] + 0.75 * Math.cos(2 * i * hx) * Math.exp(-3 * a * (k + 0.5) * ht);
                    }
                }
            }

            a_arr = new double[J + 1];
            b_arr = new double[J + 1];
            c_arr = new double[J + 1];
            d_arr = new double[J + 1];

            for(int i = 1; i < I; i++){
                m = 0;
                for(int j = 0; j < J; j++, m++){
                    if(j == 0){
                        a_arr[m] = 0;
                        b_arr[m] = 1;
                        c_arr[m] = 0;
                        d_arr[m] = Math.cos(2 * i * hx) * Math.exp(-3 * a * (k + 1) * ht);
                    }
                    else{
                        a_arr[m] = -sy;
                        b_arr[m] = 1 + 2 * sy;
                        c_arr[m] = -sy;
                        d_arr[m] = u_[i][j];
                    }
                }

                a_arr[m] = -1 / hy;
                b_arr[m] = 1 / hy;
                c_arr[m] = 0;
                d_arr[m] = 0.75 * Math.cos(2 * i * hx) * Math.exp(-3 * a * (k + 1) * ht);

                ArrayList<Double> result = Progonka(a_arr, b_arr, c_arr, d_arr);
                for(int j = 0; j < J + 1; j++){
                    u[i][j][k + 1] = result.get(j);
                }
            }
        }
        return u;
    }

    private static double U(double x, double y, double t, double a){
        return Math.cos(2 * x) * Math.cosh(y) * Math.exp(-3 * a * t);
    }

    static ArrayList<Double> Progonka(double[] a, double[] b, double[] c, double[] d)
    {
        ArrayList<Double> roots = new ArrayList<>();
        ArrayList<Double> P = new ArrayList<>();
        ArrayList<Double> Q = new ArrayList<>();

        P.add(-c[0] / b[0]);
        Q.add(d[0] / b[0]);

        for(int i = 1; i < a.length; i++)
        {
            P.add(-c[i] / (b[i] + a[i] * P.get(i - 1)));
            Q.add((d[i] - a[i] * Q.get(i - 1)) / (b[i] + a[i] * P.get(i - 1)));
        }

        Collections.reverse(P);
        Collections.reverse(Q);

        roots.add(Q.get(0));
        for (int i = 1; i < a.length; i++)
        {
            roots.add(P.get(i) * roots.get(i - 1) + Q.get(i));
        }

        Collections.reverse(roots);
        return roots;
    }
}
