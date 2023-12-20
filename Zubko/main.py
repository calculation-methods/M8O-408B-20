import math

INF = 1e10


def f(x):
    return 1.0 / (1 + x**2)


def integrate_rectangle_method(f, l, r, h):
    result = 0
    cur_x = l
    while cur_x < r:
        result += h * f((cur_x + cur_x + h) * 0.5)
        cur_x += h
    return result


def integrate_with_definite_integral(f, l, r, h=0.01, eps=1e-6):
    def f_new(t):
        return (1.0 / t**2) * f(1.0 / t)

    result = 0
    if r == INF:
        new_r = max(eps, l)
        result += integrate_rectangle_method(f_new, eps, 1.0 / new_r - eps, h)
    else:
        new_r = r
    if l == -INF:
        new_l = min(-eps, r)
        result += integrate_rectangle_method(f_new, 1.0 / new_l + eps, -eps, h)
    else:
        new_l = l
    if new_l < new_r:
        result += integrate_rectangle_method(f, new_l, new_r, h)
    return result


def integrate_lim(f, l, r, h=0.1, eps=1e-6):
    result = 0
    iters = 0
    if r == INF:
        finish = False
        cur_x = max(l, 0)
        while not finish:
            iters += 1
            new_result = result + h * f((cur_x + cur_x + h) * 0.5)
            cur_x += h
            if abs(new_result - result) < eps:
                finish = True
            result = new_result
    else:
        result += integrate_rectangle_method(f, 0, r, h)
    if l == -INF:
        finish = False
        cur_x = min(0, r)
        while not finish:
            iters += 1
            new_result = result + h * f((cur_x - h + cur_x) * 0.5)
            cur_x -= h
            if abs(new_result - result) < eps:
                finish = True
            result = new_result
    else:
        result += integrate_rectangle_method(f, l, 0, h)
    return result, iters


if __name__ == "__main__":
    a = -INF
    b = INF
    h = 0.1
    eps = 1e-9
    print("Transforming to definite integral")
    res_definite = integrate_with_definite_integral(f, a, b, h, eps)
    print("Integral =", res_definite)
    print()

    print("Limit method")
    res_limit, iters_limit = integrate_lim(f, a, b, h, eps)
    print("Integral =", res_limit)
    print("Iterations:", iters_limit)
    print()
