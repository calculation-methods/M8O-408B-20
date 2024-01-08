import numpy as np
import matplotlib.pyplot as plt
import spline as splpckg
import point as pnt
import math


class CurveFitter:
    def __init__(self, spline):
        self.m_delta = 0
        self.p = 0
        self.m_error = 0
        self.m_penalty = self.penalty(spline)

    def penalty(self, spline):
        knots = spline.get_knots()
        g = spline.get_internal_knots_num()
        self.m_penalty = 0
        for i in range(g + 1):
            self.m_penalty += 1.0 / (knots[i + 1] - knots[i])
        self.m_error = self.m_delta + self.p * self.m_penalty
        return self.m_penalty

    def delta(self, spline, points, smoothing_weight):
        self.m_delta = 0
        for point in points:
            e = point.w * (point.y - spline.get_value(point.x))
            self.m_delta += e * e
        nu = 0
        if smoothing_weight > 0:
            k = spline.get_degree()
            g = spline.get_internal_knots_num()
            coefficients = spline.get_coefficients()
            for q in range(k + 1, g + k + 1):
                e = 0
                for i in range(q - k - 1, q + 1):
                    e += coefficients[i] * spline.get_lead_derivative_difference(i, q)
                    nu += e * e
            nu *= smoothing_weight
        self.m_delta += nu
        self.m_error = self.m_delta + self.p * self.m_penalty
        return self.m_delta

    @staticmethod
    def penalty_derivative(spline, knot_id):
        knots = spline.get_knots()
        a = knots[knot_id - 1]
        b = knots[knot_id]
        c = knots[knot_id + 1]
        bma = b - a
        cmb = c - b
        return 1.0 / (cmb * cmb) - 1.0 / (bma * bma)

    def error(self, spline, points, sw):
        return self.delta(spline, points, sw) + self.p * self.penalty(spline)

    def error_gradient(self, spline, points, smoothing_weight, knot_id):
        grad_error = 0
        for point in points:
            w_sq = point.w * point.w
            diff = point.y - spline.get_value(point.x)
            grad_error -= spline.get_value_derivative_knot(point.x, knot_id + spline.get_degree()) * w_sq * diff

        if self.p > 0:
            grad_error += 0.5 * self.p * self.penalty_derivative(spline, knot_id)

        if smoothing_weight > 0:
            sm_error = 0
            k = spline.get_degree()
            g = spline.get_internal_knots_num()
            coefficients = spline.get_coefficients()
            for q in range(k + 1, g + k + 1):
                sum1 = 0
                sum2 = 0
                for i in range(q - k - 1, q + 1):
                    ci = coefficients[i]
                    lead_der_diff = spline.get_lead_derivative_difference(i, q)
                    sum1 += ci * lead_der_diff
                    sum2 += ci * spline.get_lead_der_diff_der_knot(lead_der_diff, i, q, knot_id + spline.get_degree())
                sm_error += sum1 * sum2
            grad_error += smoothing_weight * sm_error
        return 2 * grad_error

    def theta(self, spline, points, sw, alpha, direction, fixed_knots):
        g = spline.get_internal_knots_num()
        knots = [0] * (g + 2)
        knots[0] = spline.get_left_bound()
        knots[g + 1] = spline.get_right_bound()
        for i in range(g):
            knots[i + 1] = fixed_knots[i + 1] + alpha * direction[i]
        spline.set_knots(knots)

        if self.approximate(spline, points, sw):
            return self.error(spline, points, sw)
        return -1

    @staticmethod
    def norm(v):
        return sum(i * i for i in v)

    @staticmethod
    def approximate(spline, points, smoothing_weight):
        k = spline.get_degree()
        g = spline.get_internal_knots_num()
        n = len(points)
        coefficients = [0] * (g + k + 1)
        A = [[0 for i in range(g + k + 1)] for j in range(g + k + 1)]

        l = 0
        for r in range(n):
            l = spline.get_left_node_index(points[r].x, l)
            if l < 0:
                return
            b_splines = spline.b_splines(points[r].x, k)
            for i in range(k + 1):
                w_sq = points[r].w * points[r].w
                for j in range(i + 1):
                    A[i + l - k][j + l - k] += w_sq * b_splines[i] * b_splines[j]
                coefficients[i + l - k] += w_sq * points[r].y * b_splines[i]

        if smoothing_weight > 0:
            for q in range(g):
                for i in range(q, q + k + 2):
                    ai = spline.get_lead_derivative_difference(i, q + k + 1)
                    for j in range(q, i + 1):
                        A[i][j] += smoothing_weight * ai * spline.get_lead_derivative_difference(j, q + k + 1)

        for i in range(g + k + 1):
            for j in range(i):
                A[j][i] = A[i][j]

        L = np.linalg.cholesky(A)

        for i in range(g + k + 1):
            for j in range(i):
                coefficients[i] -= L[i][j] * coefficients[j]
            coefficients[i] /= L[i][i]

        for i in reversed(range(g + k + 1)):
            for j in reversed(range(i + 1, g + k + 1)):
                coefficients[i] -= L[j][i] * coefficients[j]
            coefficients[i] /= L[i][i]

        spline.set_coefficients(coefficients)

        return True

    def spec_dimensional_minimization(self, spline, points, sw, direction, error_derivative, fixed_knots):
        g = spline.get_internal_knots_num()
        knots = spline.get_knots()
        alpha_max = math.inf
        a = spline.get_left_bound()
        b = spline.get_right_bound()

        if direction[0] < 0:
            alpha_max = (a - knots[1]) / direction[0]
        for i in range(g - 1):
            if direction[i] > direction[i + 1]:
                alpha_max = min(alpha_max, (knots[i + 2] - knots[i + 1]) / (direction[i] - direction[i + 1]))
        if direction[g - 1] > 0:
            alpha_max = min(alpha_max, (b - knots[g]) / direction[g - 1])

        theta0 = self.m_error
        theta0_der = 0
        for dir, error_der in zip(direction, error_derivative):
            theta0_der += dir * error_der

        alpha0 = 0
        alpha2 = alpha_max / (1 - theta0 / alpha_max / theta0_der)
        alpha1 = 0.5 * alpha2
        q0 = self.m_delta
        r0 = self.m_penalty
        theta1 = self.theta(spline, points, sw, alpha1, direction, fixed_knots)
        if theta1 < 0:
            return False
        q1 = self.m_delta
        r1 = self.m_penalty

        iteration = 0
        max_num_of_iterations = 10
        while theta1 >= theta0 and iteration < max_num_of_iterations:
            alpha_tilde = -0.5 * theta0_der * alpha1 * alpha1 / (theta1 - theta0 - theta0_der * alpha1)
            alpha1 = max(0.1 * alpha1, alpha_tilde)
            theta1 = self.theta(spline, points, sw, alpha1, direction, fixed_knots)
            if theta1 < 0:
                return False
            q1 = self.m_delta
            r1 = self.m_penalty
            iteration += 1

        if iteration > 0:
            if theta1 > theta0:
                self.theta(spline, points, sw, alpha0, direction, fixed_knots)
            return True

        theta2 = self.theta(spline, points, sw, alpha2, direction, fixed_knots)
        if theta2 < 0:
            return False
        q2 = self.m_delta
        r2 = self.m_penalty

        while theta2 < theta1:
            alpha0 = alpha1
            q0 = q1
            r0 = r1
            alpha1 = alpha2
            theta1 = theta2
            q1 = q2
            r1 = r2

            alpha2 = min(2 * alpha1, 0.5 * (alpha_max + alpha1))
            theta2 = self.theta(spline, points, sw, alpha2, direction, fixed_knots)
            if theta2 < 0:
                return False

            q2 = self.m_delta
            r2 = self.m_penalty

        a0 = q0
        diff1 = alpha1 - alpha0
        diff2 = alpha2 - alpha0
        a2 = (q1 - q0) / diff1
        a2 -= (q2 - q0) / diff2
        a2 /= alpha1 - alpha2
        a1 = (q1 - a0) / diff1
        a1 -= a2 * diff1

        fraction = diff1 / diff2
        numerator = r1 - r0 - fraction * (r2 - r0)
        temp = math.log((alpha_max - alpha1) / (alpha_max - alpha0))
        denominator = temp - fraction * math.log((alpha_max - alpha2) / (alpha_max - alpha0))
        b2 = numerator / denominator
        b1 = (r1 - r0 - b2 * temp) / diff1

        a = -2 * a2
        b = -a * (alpha_max + alpha0) - self.p * b1 - a1
        c = (self.p * b1 + a1 + a * alpha0) * alpha_max - self.p * b2

        root1 = -0.5 * (b + math.sqrt(b * b - 4 * a * c)) / a
        root2 = -b / a - root1

        alpha_res = 0
        if 0 < root1 < alpha_max:
            alpha_res = root1
        elif 0 < root2 < alpha_max:
            alpha_res = root2

        theta_res = self.theta(spline, points, sw, alpha_res, direction, fixed_knots)

        if theta_res < 0:
            return False
        return True

    @staticmethod
    def initiate_grid(spline, points):
        k = spline.get_degree()
        g = spline.get_internal_knots_num()
        knots = [0] * (g + 2)
        knots[0] = spline.get_left_bound()
        knots[g + 1] = spline.get_right_bound()
        n = len(points)

        unique_size = 0
        index = 0

        while index < n and points[index].x < knots[g + 1]:
            if index != 0 and points[index].x != points[index - 1].x:
                unique_size += 1
            index += 1

        if unique_size <= 0:
            return False

        if unique_size < g + k + 1:
            return False

        points_per_knot = unique_size / (g + 1)
        knot_index = 1
        i = 1
        counter = 0

        while knot_index < g + 1:
            while counter < knot_index * points_per_knot or points[i].x == points[i - 1].x:
                if points[i].x != points[i - 1].x:
                    counter += 1
                i += 1
            knots[knot_index] = 0.5 * (points[i].x + points[i - 1].x)
            knot_index += 1

        spline.set_knots(knots)
        return True

    def approximate_with_optimal_grid(self, spline, points, smooth_weight, eps1, eps2):
        if not self.initiate_grid(spline, points):
            return False

        g = spline.get_internal_knots_num()

        if not self.approximate(spline, points, smooth_weight):
            return False

        direction = [0] * g
        error_derivative = [0] * g

        self.delta(spline, points, smooth_weight)
        self.p = eps1 * self.m_delta * (spline.get_right_bound() - spline.get_left_bound()) / (g + 1) / (g + 1)
        self.m_error = self.m_delta + self.p * self.penalty(spline)

        for i in range(g):
            error_derivative[i] = self.error_gradient(spline, points, smooth_weight, i + 1)
            direction[i] = -error_derivative[i]

        old_norm = self.norm(direction)
        criteria1 = eps1 + eps2
        criteria2 = criteria1
        max_num_of_iter = 1000

        iteration = 0
        eps2_sq = eps2 * eps2

        while (criteria1 >= eps1 or criteria2 >= eps2_sq) and iteration < max_num_of_iter:
            fixed_knots = spline.get_knots().copy()

            old_error = self.m_error

            if not self.spec_dimensional_minimization(spline, points, smooth_weight, direction, error_derivative, fixed_knots):
                spline.set_knots(fixed_knots)
                return self.approximate(spline, points, smooth_weight)

            for i in range(g):
                error_derivative[i] = self.error_gradient(spline, points, smooth_weight, i + 1)

            new_norm = self.norm(error_derivative)

            if iteration % g == 0:
                for i in range(g):
                    direction[i] = -error_derivative[i]
            else:
                temp = new_norm / old_norm
                for i in range(g):
                    direction[i] *= temp
                    direction[i] -= error_derivative[i]

            numerator = 0
            denominator = 0

            knots = spline.get_knots()
            for i in range(1, len(knots) - 1):
                temp = knots[i] - fixed_knots[i]
                numerator += temp * temp
                denominator += fixed_knots[i] * fixed_knots[i]

            criteria1 = math.fabs(old_error - self.m_error) / old_error
            criteria2 = numerator / denominator

            old_norm = new_norm

        return True


def get_x(length):
    x = np.zeros(length)
    for i in range(length):
        x[i] = i + 1.0 / (i + 1) * np.random.rand()
    x = np.concatenate((x, x), 0)
    x = np.sort(x)

    return x


def get_data_1():
    x = get_x(60)
    n = len(x)
    y = [0] * n
    w = [1] * n
    for i in range(0, n, 2):
        y[i] = np.cos(0.2 * x[i])
        err = np.random.rand() * 15 / (x[i] + 10)
        y[i + 1] = y[i] + err
        y[i] -= err
        w[i] = 1.0 / math.fabs(x[i] - x[i - 1])
        w[i + 1] = w[i]

    return pnt.Points(x, y, w)


def get_data_2():
    x = get_x(60)
    n = len(x)
    y = [0] * n
    w = [1] * n
    for i in range(0, n, 2):
        y[i] = np.cos(0.2 * x[i])
        y[i] += 0.4 * np.cos(0.5 * x[i])
        err = np.random.rand() / 3
        y[i + 1] = y[i] + err
        y[i] -= err
        w[i] = 1.0 / math.fabs(x[i] - x[i - 1])
        w[i + 1] = w[i]

    return pnt.Points(x, y, w)


def get_data_3():
    x = get_x(60)
    n = len(x)
    y = [0] * n
    w = [1] * n
    for i in range(0, n, 2):
        y[i] = np.cos(0.005 * x[i] * x[i])
        err = np.random.rand() * 15 / (x[i] + 10)
        y[i + 1] = y[i] + err
        y[i] -= err
        w[i] = 1.0 / math.fabs(x[i] - x[i - 1])
        w[i + 1] = w[i]

    return pnt.Points(x, y, w)


def main():
    points = get_data_2()

    x_curve = np.linspace(points.x[0], points.x[-1], 1000)
    y_curve_opt = [None] * 1000

    knots = [points.x[0], 1, 2, 3, 4, points.x[len(points) - 1]]
    coefficients = [None] * (len(points) + 2)
    knots.append(knots[len(knots) - 2] + 1)
    knots.sort()
    s = splpckg.Spline(coefficients, knots, 3)
    c = CurveFitter(s)
    c.initiate_grid(s, points)
    q = 1e-9

    c.approximate_with_optimal_grid(s, points, q, 1e-3, 1e-3)

    index = 0
    for point in x_curve:
        y_curve_opt[index] = s.get_value(point)
        index += 1

    index = 0
    knots2 = s.get_knots()
    knots_y = [0] * len(knots2)
    for point in knots2:
        knots_y[index] = s.get_value(point)
        index += 1

    plt.figure(figsize=(10, 6))
    plt.plot(x_curve, y_curve_opt, points.x, points.y, 'o', knots2, knots_y, 's')
    plt.legend(["Сплайн", "Данные", "Оптимальные узлы"])
    plt.xlim([points.x[0] - 1, points.x[-1] + 1])

    plt.savefig('sample.png')

main()