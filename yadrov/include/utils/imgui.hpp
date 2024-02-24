#pragma once

template <class... Args>
std::unique_ptr<char[]> format(const char* str, Args... args) {
    int size = std::snprintf(nullptr, 0, str, args...) + 1;
    if (size <= 0) { 
    	throw std::runtime_error("Error during formatting." );
    }
    std::unique_ptr<char[]> buf(new char[size]);
    std::snprintf(buf.get(), size, str, args...);
    return buf;
}

struct AxisSharedLimits {
	double x_min, x_max, y_min, y_max;

	void assign(double x_min_, double x_max_, double y_min_, double y_max_, double padding = 0) {
		assert(x_min_ <= x_max_ && y_min_ <= y_max_);
		double dx = x_max_ - x_min_, dy = y_max_ - y_min_;
		x_min = x_min_ - dx * padding;
		x_max = x_max_ + dx * padding;
		y_min = y_min_ - dy * padding;
		y_max = y_max_ + dy * padding;
	}

	void apply() {
		ImPlot::SetupAxisLinks(ImAxis_X1, &x_min, &x_max);
		ImPlot::SetupAxisLinks(ImAxis_Y1, &y_min, &y_max);
	}
};