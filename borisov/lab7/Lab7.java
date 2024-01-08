import java.util.ArrayList;

public class Lab7 {
    static int count = 0;
    public static void main(String[] args) {
        double h = 0.1;
        int flag = 0;
        double[][] u = SolveEqLaplas(h, flag);
        int s = 9;

        double max = 0;
        double delta = 0;
        double temp;
        for (int i = 0; i < u[0].length; i++) {
            temp = Math.abs(u[s][i] - U(h * i, h * s));
            delta += temp * temp;
            if(temp > max){
                max = temp;
            }
        }
        System.out.println("Метод простых итераций:");
        System.out.println(delta);
        System.out.println("Число шагов: " + count);

        flag = 1;
        u = SolveEqLaplas(h, flag);

        max = 0;
        delta = 0;
        for (int i = 0; i < u[0].length; i++) {
            temp = Math.abs(u[s][i] - U(h * i, h * s));
            delta += temp * temp;
            if(temp > max){
                max = temp;
            }
        }
        System.out.println("\nМетод Зейделя:");
        System.out.println(delta);
        System.out.println("Число шагов: " + count);

        flag = 2;
        u = SolveEqLaplas(h, flag);

        max = 0;
        delta = 0;
        for (int i = 0; i < u[0].length; i++) {
            temp = Math.abs(u[s][i] - U(h * i, h * s));
            delta += temp * temp;
            if(temp > max){
                max = temp;
            }
        }
        System.out.println("\nМетод простых итераций с верхней релаксацией:");
        System.out.println(delta);
        System.out.println("Число шагов: " + count);
    }

    private static double[][] SolveEqLaplas(double h, int flag){
        int N = (int)(Math.PI / h);
        int M = (int)(1 / h);
        double[] column = new double[(N + 1) * (M + 1)];
        double[][] A = new double[(N + 1) * (M + 1)][(N + 1) * (M + 1)];

        int k = 1;
        int l = 1;
        int p = 0;
        for(int i = 0; i < A.length; i++){
            for(int j = 0; j < A[0].length; j++){
                if(i >= 0 && i < N + 1){
                    if(j == i) {
                        A[i][j] = 1;
                        column[i] = Math.sin(i * h);
                    }
                }
                else if(i == k * (N + 1) && k < M){
                    if(j == i) {
                        A[i][j] = 1;
                        column[i] = -h * Math.exp(k * h);
                    }
                    else if(j - 1 == i) A[i][j] = -1;
                    else if(j == A[0].length - 1) k++;
                }
                else if(i == l * (N + 1) + N && l < M){
                    if(j == i) {
                        A[i][j] = 1;
                        column[i] = -h * Math.exp(l * h);
                    }
                    else if(j + 1 == i) A[i][j] = -1;
                    else if(j == A[0].length - 1) l++;
                }
                else if(i >= k * (N + 1)){
                    if(j == i) {
                        A[i][j] = 1;
                        column[i] = Math.E * Math.sin(p * h);
                        p++;
                    }
                }
                else{
                    if(j == i) {
                        A[i][j] = -4;
                        column[i] = 0;
                    }
                    else if(j - 1 == i) A[i][j] = 1;
                    else if(j + 1 == i) A[i][j] = 1;
                    else if(j + (N + 1) == i) A[i][j] = 1;
                    else if(j - (N + 1) == i) A[i][j] = 1;
                }
            }
        }

        double[][] u_ = new double[M + 1][N + 1];
        if(flag == 0) {
            double[] u = Iteration(A, column);
            int n = 0;
            for (int i = 0; i < u_.length; i++) {
                for (int j = 0; j < u_[0].length; j++) {
                    u_[i][j] = u[n];
                    n++;
                }
            }
        }
        else if(flag == 1){
            double[] u = Zeidel(A, column);
            int n = 0;
            for (int i = 0; i < u_.length; i++) {
                for (int j = 0; j < u_[0].length; j++) {
                    u_[i][j] = u[n];
                    n++;
                }
            }
        }
        else if(flag == 2){
            double[] u = IterationWithRelax(A, column);
            int n = 0;
            for (int i = 0; i < u_.length; i++) {
                for (int j = 0; j < u_[0].length; j++) {
                    u_[i][j] = u[n];
                    n++;
                }
            }
        }
        return u_;
    }


    private static double U(double x, double y){
        return Math.sin(x) * Math.exp(y);
    }
    private static double[] Iteration(double[][] matrix, double[] column)
    {
        double[] result = new double[column.length];

        double[][] alpha = new double[matrix.length][matrix[0].length];
        double[] betta = new double[column.length];
        double[] x_cur = new double[column.length];
        double[] x_prev = new double[column.length];
        boolean norm = true;
        double sum = 0;
        double max_sum = 0;
        int count_iteration = 0;

        for (int i = 0; i < matrix.length; i++)
        {
            sum = 0;
            for (int j = 0; j < matrix[0].length + 1; j++)
            {
                if (j != i && j != matrix[0].length)
                {
                    alpha[i][j] = -matrix[i][j] / matrix[i][i];
                    sum += Math.abs(alpha[i][j]);
                }
                else if (j == i) alpha[i][j] = 0;

                if (j == matrix[0].length)
                {
                    betta[i] = column[i] / matrix[i][i];
                    x_prev[i] = betta[i];
                }
            }

            if (sum > max_sum) max_sum = sum;
        }

        if (max_sum > 1) norm = false;

        if(norm)
        {
            double epsilon = 0.0001;
            double epsilon_i = 1;

            while (epsilon_i > epsilon)
            {
                epsilon_i = 0;
                x_cur = SumVectors(betta, MultyMatrVector(alpha, x_prev));

                for (int i = 0; i < column.length; i++)
                    epsilon_i += Math.pow(x_prev[i] - x_cur[i], 2);
                epsilon_i = Math.sqrt(epsilon_i);

                x_prev = x_cur;
                count_iteration++;
            }
        }
        count = count_iteration;
        result = x_cur;
        return result;
    }

    static double[] MultyMatrVector(double[][] matrix, double[] column)
    {
        double[] result = new double[column.length];
        for(int i = 0; i < matrix.length; i++)
        {
            for(int j = 0; j < matrix[0].length; j++)
            {
                result[i] += matrix[i][j] * column[j];
            }
        }
        return result;
    }

    static double[] SumVectors(double[] a, double[] b)
    {
        double[] result = new double[a.length];

        for (int i = 0; i < a.length; i++)
        {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    static double[] MultyNumberVector(double[] a, double lambda){
        double[] res = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            res[i] = lambda * a[i];
        }
        return res;
    }

    static double[] Zeidel(double[][] matrix, double[] column)
    {
        double[] result = new double[column.length];

        double[][] alpha = new double[matrix.length][matrix[0].length];
        double[] betta = new double[column.length];
        double[] x_prev = new double[column.length];
        double[] x_cur = new double[column.length];
        boolean norm = true;
        double sum = 0;
        double max_sum = 0;
        int count_iteration = 0;

        for (int i = 0; i < matrix.length; i++)
        {
            sum = 0;
            for (int j = 0; j < matrix[0].length + 1; j++)
            {
                if (j != i && j != matrix[0].length)
                {
                    alpha[i][j] = -matrix[i][j] / matrix[i][i];
                    sum += Math.abs(alpha[i][j]);
                }
                else if (j == i) alpha[i][j] = 0;

                if (j == matrix[0].length)
                {
                    betta[i] = column[i] / matrix[i][i];
                    x_prev[i] = betta[i];
                }
            }

            if (sum > max_sum) max_sum = sum;
        }
        if (max_sum > 1) norm = false;

        double[] vctr = new double[column.length];

        if(norm)
        {
            double epsilon = 0.0001;
            double epsilon_i = 1;

            ArrayList<double[]> str = StrOfMatr(alpha);
            while (epsilon_i > epsilon)
            {
                epsilon_i = 0;

                for (int i = 0; i < vctr.length; i++)
                {
                    vctr[i] = x_prev[i];
                }

                for (int i = 0; i < x_cur.length; i++)
                {
                    x_cur[i] = betta[i] + MultyStrVector(str.get(i), x_prev);
                    x_prev[i] = x_cur[i];
                }

                for (int i = 0; i < column.length; i++)
                    epsilon_i += Math.pow(vctr[i] - x_cur[i], 2);
                epsilon_i = Math.sqrt(epsilon_i);
                count_iteration++;
            }

        }
        count = count_iteration;
        result = x_cur;
        return result;
    }

    static ArrayList<double[]> StrOfMatr(double[][] matrix)
    {
        ArrayList<double[]> str = new ArrayList<>();
        //double[] mas = new double[matrix.GetLength(0)];

        for(int i = 0; i < matrix.length; i++)
        {
            double[] mas = new double[matrix.length];
            for (int j = 0; j < matrix[0].length; j++)
            {
                mas[j] = matrix[i][j];
            }

            str.add(mas);
        }
        return str;
    }

    static double MultyStrVector(double[] str, double[] vctr)
    {
        double result = 0;

        for(int i = 0; i < str.length; i++)
        {
            result += str[i] * vctr[i];
        }
        return result;
    }

    static double[] IterationWithRelax(double[][] matrix, double[] column){
        double[] result = new double[column.length];

        double[][] alpha = new double[matrix.length][matrix[0].length];
        double[] betta = new double[column.length];
        double[] x_cur = new double[column.length];
        double[] x_prev = new double[column.length];
        double[] x_predict = new double[column.length];
        boolean norm = true;
        double sum = 0;
        double max_sum = 0;
        int count_iteration = 0;
        double w = 1.01;

        for (int i = 0; i < matrix.length; i++)
        {
            sum = 0;
            for (int j = 0; j < matrix[0].length + 1; j++)
            {
                if (j != i && j != matrix[0].length)
                {
                    alpha[i][j] = -matrix[i][j] / matrix[i][i];
                    sum += Math.abs(alpha[i][j]);
                }
                else if (j == i) alpha[i][j] = 0;

                if (j == matrix[0].length)
                {
                    betta[i] = column[i] / matrix[i][i];
                    x_prev[i] = betta[i];
                }
            }

            if (sum > max_sum) max_sum = sum;
        }

        if (max_sum > 1) norm = false;

        if(norm)
        {
            double epsilon = 0.0001;
            double epsilon_i = 1;

            while (epsilon_i > epsilon)
            {
                epsilon_i = 0;
                x_predict = SumVectors(betta, MultyMatrVector(alpha, x_prev));
                x_cur = SumVectors(MultyNumberVector(x_predict, w), MultyNumberVector(x_prev, 1 - w));

                for (int i = 0; i < column.length; i++)
                    epsilon_i += Math.pow(x_prev[i] - x_cur[i], 2);
                epsilon_i = Math.sqrt(epsilon_i);

                x_prev = x_cur;
                count_iteration++;
            }
        }
        count = count_iteration;
        result = x_cur;
        return result;
    }
}
