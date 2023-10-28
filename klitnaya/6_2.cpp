#include <iostream>
#include <cmath>
#include <vector>


//аналитическое решение
double analitic(double a, double x, double t) {
    return std::sin(x - a * t);
}

//начальные условия
double u_x_0(double x) {
    return std::sin(x);
}

double u_x_0t(double x, double a) {
    return -a * std::cos(x);
}


//граничные условия
double u_l(double a, double t) {
    return -std::sin(a * t);
}

double u_r(double a, double t) {
    return std::sin(a * t);
}


//аналитическое решение на сетке
std::vector<std::vector<double>> analitic_on_net(double t_, int lc, double h, double a) {
    int x_c = static_cast<int>(M_PI / h) + 1;
    std::vector<std::vector<double>> group(lc, std::vector<double>(x_c, 0.0));
    for (int i = 0; i < lc; i++) {
        for (int j = 0; j < x_c; j++) {
            group[i][j] = analitic(a, j * h, i * t_);
        }
    }
    return group;
}

//явный метод
 std::vector<std::vector<double>> explicit_solution(double t_, int lc, double h, int order, double a) {
    int x_count = static_cast<int>(M_PI / h) + 1;
    std::vector<std::vector<double>> group(lc, std::vector<double>(x_count, 0.0));
    for (int i = 0; i < lc; i++) {
        if (i == 0) {
            for (int j = 0; j < x_count; j++) {
                group[i][j] = u_x_0(j * h);
            }
        }
        else if (i == 1) {
            if (order == 1) {
                for (int j = 0; j < x_count; j++) {
                    group[i][j] = group[i-1][j] + u_x_0t(j * h, a) * t_;
                }
            }
            else if (order == 2) {
                for (int j = 0; j < x_count; j++) {
                    group[i][j] = group[i-1][j] + u_x_0t(j * h, a) * t_ + (a * (-u_x_0(j * h))) * (t_ * t_) / 2;
                }
            }
        }
        else {
            group[i][0] = u_l(i * t_, a);
            for (int j = 1; j < x_count - 1; j++) {
                group[i][j] = (t_ * t_) / (h * h) * (group[i-1][j+1] - 2 * group[i-1][j] + group[i-1][j-1]) + 2 * group[i-1][j] - group[i-2][j];
            }
            group[i][x_count-1] = u_r(i * t_, a);
        }
    }
    return group;
}
 
//неявный метод
std::vector<std::vector<double>> implicit_solve(double t_, int lc, double h, int order, double a) {
    int x_count = static_cast<int>(M_PI / h) + 1;
    std::vector<std::vector<double>> group(lc, std::vector<double>(x_count, 0.0));
    std::vector<double> alpha(x_count - 1, 0.0);
    std::vector<double> beta(x_count - 1, 0.0);
    double aa = -a * (t_ * t_) / (h * h);
    double bb = 1 + 2 * a * (t_ * t_) / (h * h);
    double cc = aa;
    for (int i = 0; i < lc; i++) {
        if (i == 0) {
            for (int j = 0; j < x_count; j++) {
                group[i][j] = u_x_0(j * h);
            }
        }
        else if (i == 1) {
            if (order == 1) {
                for (int j = 0; j < x_count; j++) {
                    group[i][j] = group[i-1][j] + u_x_0t(j * h, a) * t_;
                }
            }
            else if (order == 2) {
                for (int j = 0; j < x_count; j++) {
                    group[i][j] = group[i-1][j] + u_x_0t(j * h, a) * t_ + (a * (-u_x_0(j * h))) * (t_ * t_) / 2;
                }
            }
        }
        else {
            for (int j = 0; j < x_count; j++) {
                group[i][j] = 0.0;
            }
            beta[0] = u_l(i*t_, a);
            for (int j = 1; j < x_count - 1; j++) {
                alpha[j] = -aa / (bb + cc * alpha[j-1]);
                beta[j] = (2 * group[i-1][j] - group[i-2][j] - cc * beta[j-1]) / (bb + cc * alpha[j-1]);
            }
            group[i][x_count - 1] = u_r(i * t_, a);
            for (int j = x_count - 2; j >= 0; j--) {
                group[i][j] = group[i][j+1] * alpha[j] + beta[j];
            }
        }
    }
    return group;
}

//поиск ошибки
double error_ (std::vector<double>& explicit_result, std::vector<double>& analitic_resul){
    double er = 0.0;
    for (int i = 0; i < explicit_result.size(); i++){
        er += (explicit_result[i] - analitic_resul[i]) * (explicit_result[i] - analitic_resul[i]);
    }
    er = er / explicit_result.size();
    return er;
  }




int main() {
    // Пример использования функции
    double t_ = 0.01;
    int lc = 5;
    double h = 0.1;
    int order = 2;
    double a = 0.5;
    std:: cout << "Введите коэффициент а: ";
    std:: cin >>a;
    std:: cout << "Введите шаг по временной шкале: ";
    std:: cin >> t_;
    std:: cout << "Введите шаг по пространственной шкале: ";
    std:: cin >> h;
    if (h < sqrt(a) * t_ ){
        char tmp;
        std:: cout << "Вы ввели неверный пространственный шаг, поэтому явный метод будет неустойчив. Если хотите поменять параметр нажмите y иначе n ";
        std:: cin >>tmp;
        if (tmp == 'y'){
            std::cin >> h;
        }
    }
    std:: cout <<"Введите порядок аппроксимации: ";
    std:: cin >> order;
    //std::cin >> t_ >> lc >> h >> order >> a; 
    
    
    std::vector<std::vector<double>> explicit_sol = explicit_solution(t_, lc, h, order, a);
    std::vector<std::vector<double>> implicit_sol = implicit_solve(t_, lc, h, order, a);
    std::vector<std::vector<double>> analitic_sol = analitic_on_net(t_, lc, h, a);
    
    
    std::cout << "Analitic_solution Explicit_solution Implicit_solution: "<<std::endl;
    for (int i = 0; i < explicit_sol.size() ; i++) {
        for (int j = 0; j < explicit_sol[0].size() ; j++) {
            std:: cout <<analitic_sol[i][j] << " " <<explicit_sol[i][j]<< " " << implicit_sol[i][j]<<std::endl;
        }
    }
    std::cout <<std::endl;
    std:: cout << "error for:  Explicit_solution Implicit_solution" << std::endl;
    for (int i = 0; i < explicit_sol.size() ; i++){
        double t1 = error_(explicit_sol[i], analitic_sol[i]);
        double t2 = error_(implicit_sol[i], analitic_sol[i]);
        std::cout << t1 << " "<< t2<<std::endl;
    }
    return 0;
}


