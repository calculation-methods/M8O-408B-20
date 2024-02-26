#include "imgui.h"
#include "imgui_internal.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"
#include "implot.h"
#include <GLFW/glfw3.h>

#include <iostream>
#include <algorithm>
#include <memory>
#include <cstdio>

#include "../../include/partial_differential/parabolic_pde.hpp"
#include "presets.hpp"

using namespace ParabolicPDE;

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

void ShowSolutionWindow(const char* name, const std::vector<double>& x, const std::vector<double>& t, const grid_t<double>& u_numeric, 
                        const std::vector<double>& error, const std::vector<int>& timestamps, const grid_t<double>& u_correct, bool show_correct_solution, 
                        AxisSharedLimits& solution_limits, AxisSharedLimits& error_limits, ImGuiWindowFlags flags = 0) {

	ImGui::Begin(name, nullptr, flags);
	if (ImPlot::BeginPlot("Solutions")) {
		ImPlot::PushColormap("Paired");
		ImPlot::SetupAxes("x", "u(x)");
		solution_limits.apply();
		if (!u_numeric.empty() && !timestamps.empty()) {
			for (size_t i = 0; i < timestamps.size(); ++i) {
				if (show_correct_solution) {
					ImPlot::PushColormap("Deep");
					ImPlot::SetNextLineStyle(ImPlot::GetColormapColor(8));
					ImPlot::PlotLine("", x.data(), u_correct[timestamps[i]].data(), x.size());
					ImPlot::PopColormap();
				}
				ImPlot::SetNextLineStyle(ImPlot::GetColormapColor(i));
				ImPlot::PlotLine(format("t = %.3f##%d", t[timestamps[i]], i).get(), x.data(), u_numeric[timestamps[i]].data(), x.size());
			}
		}
		ImPlot::PopColormap();
		ImPlot::EndPlot();
	}

	if (ImPlot::BeginPlot("Error")) {
		ImPlot::PushColormap("Deep");
		ImPlot::SetupAxes("t", "err");
		error_limits.apply();
		ImPlot::PlotLine("", t.data(), error.data(), t.size());
		ImPlot::PopColormap();
		ImPlot::EndPlot();
	}
	ImGui::End();
}

static void glfw_error_callback(int error, const char* description) {
	std::cerr << "Glfw Error "<< error << ": " << description << '\n';
}

int main() {
	// Setup window
	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit())
		return 1;

	// Decide GL+GLSL versions
#if defined(__APPLE__)
	// GL 3.2 + GLSL 150
	const char* glsl_version = "#version 150";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
#else
	// GL 3.0 + GLSL 130
	const char* glsl_version = "#version 130";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
#endif

	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	GLFWwindow* window = glfwCreateWindow(1550, 900, "lab5", NULL, NULL);
	if (window == NULL)
		return 1;
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImPlot::CreateContext();
	ImGuiIO& io = ImGui::GetIO();
	io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;

	ImGui::StyleColorsLight();

	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(glsl_version);

	ImVec4 clear_color = ImVec4(0, 0, 0, 1);

	std::vector<Preset> presets;
	GeneratePresets(presets);

	int idx = 1;
	PDE<double>* pde = &presets[idx].pde;
	int h_count = presets[idx].h_count;
	double t_end = presets[idx].t_end;

	ApproxType approx_type = ApproxType::Linear;
	double sigma = 0.5, theta = 0.5;
	std::vector<double> t, x, error_impl, error_expl, error_comb;
	std::vector<int> timestamps{0};
	grid_t<double> u_impl, u_expl, u_comb, u_correct;
	bool show_correct_solution = true, always_reset_axis = true;
	double max_error = 1, max_u = 1, min_u = 0;
	int tau_count = CourantCondition(h_count, sigma, t_end, pde->end, pde->a);
	AxisSharedLimits solution_limits{pde->start, pde->end, 0.0, 1.0}, error_limits{0.0, 1.0, 0.0, 1.0};

	ImGuiWindowFlags dockspace_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
	                                   ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoDocking,
                     window_flags = ImGuiWindowFlags_NoMove;
	bool first_loop = true;

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		const ImGuiViewport* viewport = ImGui::GetMainViewport();
		ImGui::SetNextWindowPos(viewport->Pos);
		ImGui::SetNextWindowSize(viewport->Size);
		ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
		ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);

		{
			ImGui::Begin("Dockspace", nullptr, dockspace_flags);
			ImGui::PopStyleVar(2);
			ImGuiID dockspace_id = ImGui::GetID("MainWindowGroup");
			ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f));
			if (first_loop) {
				first_loop = false;
				ImGui::DockBuilderRemoveNode(dockspace_id);
				ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
				ImGui::DockBuilderSetNodeSize(dockspace_id, viewport->Size);
	
				ImGuiID dock1 = ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 0.2, nullptr, &dockspace_id);
				ImGuiID dock2 = ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 1.0 / 3, nullptr, &dockspace_id);
				ImGuiID dock3 = ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 0.5, nullptr, &dockspace_id);
				ImGuiID dock4 = dockspace_id;
	
				ImGui::DockBuilderDockWindow("Parameters", dock1);
				ImGui::DockBuilderDockWindow("Presets", dock1);
				ImGui::DockBuilderDockWindow("Explicit", dock2);
				ImGui::DockBuilderDockWindow("Implicit", dock3);
				ImGui::DockBuilderDockWindow("Combined", dock4);
				ImGui::DockBuilderFinish(dockspace_id);
			}
			ImGui::End();
		}

		{
			ImGui::Begin("Presets");
			bool changed = false;
			for (size_t i = 0; i < presets.size(); ++i) {
				if (ImGui::Selectable(presets[i].name.c_str(), idx == static_cast<int>(i))) {
					idx = i;
					changed = true;
				}
			}
			ImGui::End();

			if (changed) {
				pde = &presets[idx].pde;
				h_count = presets[idx].h_count;
				t_end = presets[idx].t_end;
			}
		}

		{
			ImGui::Begin("Parameters", nullptr, window_flags);

			// TODO: Заменить все на Drag*, Slider* (imgui_demo.cpp:2075)

			ImGui::Text("Equation");
			if (ImGui::InputDouble("a", &pde->a, 1, 1, "%.2f")) {
				tau_count = CourantCondition(h_count, sigma, t_end, pde->end, pde->a);
			}
			ImGui::InputDouble("b", &pde->b, 1, 1, "%.2f");
			ImGui::InputDouble("c", &pde->c, 1, 1, "%.2f");
			ImGui::Text(format("f(x, t) = %s", presets[idx].f.c_str()).get());
			ImGui::Separator();

			ImGui::Text("Boundary conditions");
			ImGui::InputDouble("start", &pde->start, 0.1, 0.1, "%.2f");
			ImGui::InputDouble("end", &pde->end, 0.1, 0.1, "%.2f");
			ImGui::Text(format("u(x, 0) = %s", presets[idx].psi.c_str()).get());
			ImGui::Text(format("alpha1 = %.2lf, beta1 = %.2lf", pde->alpha1, pde->beta1).get());
			ImGui::Text(format("gamma1(t) = %s", presets[idx].gamma1.c_str()).get());
			ImGui::Text(format("alpha2 = %.2lf, beta2 = %.2lf", pde->alpha2, pde->beta2).get());
			ImGui::Text(format("gamma2(t) = %s", presets[idx].gamma2.c_str()).get());
			ImGui::Separator();

			ImGui::Text("Methods parameters");
			ImGui::InputDouble("End time", &t_end, 0.1, 0.1, "%.2f");
			if (ImGui::InputInt("X intervals", &h_count)) {
				tau_count = CourantCondition(h_count, sigma, t_end, pde->end, pde->a);
			}
			if (ImGui::InputDouble("Sigma", &sigma, 0.1, 0.1, "%.2f")) {
				tau_count = CourantCondition(h_count, sigma, t_end, pde->end, pde->a);
			}
			ImGui::Text(format("T intervals: %d", tau_count).get());
			ImGui::InputDouble("Theta", &theta, 0.1, 0.1, "%.2f");
			ImGui::Separator();

			ImGui::Text("Approximation");
			ImGui::RadioButton("Linear", reinterpret_cast<int*>(&approx_type), static_cast<int>(ApproxType::Linear));
			ImGui::RadioButton("Quadratic", reinterpret_cast<int*>(&approx_type), static_cast<int>(ApproxType::Quadratic));
			ImGui::RadioButton("Taylor", reinterpret_cast<int*>(&approx_type), static_cast<int>(ApproxType::Taylor));
			ImGui::Separator();

			ImGui::Text(format("u(x, t) = %s", presets[idx].solution.c_str()).get());
			ImGui::Checkbox("Show analitical solution", &show_correct_solution);
			if (ImGui::Button("Reset axis limits")) {
				solution_limits.assign(pde->start, pde->end, min_u, max_u, 0.05);
				error_limits.assign(0, t_end, 0, max_error, 0.05);
			}
			ImGui::SameLine();
			ImGui::Checkbox("Always reset", &always_reset_axis);
			if (ImGui::Button("Solve")) {
				std::tie(x, t, u_expl) = ExplicitSolver(*pde, t_end, h_count, sigma, approx_type);
				u_impl = ImplicitSolver(*pde, x, t, t_end, approx_type);
				u_comb = CombinedSolver(*pde, x, t, t_end, approx_type, theta);
				
				u_correct.assign(t.size(), std::vector<double>(x.size()));
				error_expl.assign(t.size(), 0);
				error_impl.assign(t.size(), 0);
				error_comb.assign(t.size(), 0);
				max_error = 0;
				max_u = u_impl[0][0];
				min_u = u_impl[0][0];

				for (size_t k = 0; k < t.size(); ++k) {
					double max_error_expl = 0, max_error_impl = 0, max_error_comb = 0;
					for (size_t i = 0; i < x.size(); ++i) {
						u_correct[k][i] = pde->solution(x[i], t[k]);
						max_error_expl = std::max(max_error_expl, std::abs(pde->solution(x[i], t[k]) - u_expl[k][i]));
						max_error_impl = std::max(max_error_impl, std::abs(pde->solution(x[i], t[k]) - u_impl[k][i]));
						max_error_comb = std::max(max_error_comb, std::abs(pde->solution(x[i], t[k]) - u_comb[k][i]));
						
						max_u = std::max({max_u, u_expl[k][i], u_impl[k][i], u_comb[k][i]});
						min_u = std::min({min_u, u_expl[k][i], u_impl[k][i], u_comb[k][i]});
					}
					error_expl[k] = max_error_expl;
					error_impl[k] = max_error_impl;
					error_comb[k] = max_error_comb;
					max_error = std::max({max_error, max_error_impl, max_error_expl, max_error_comb});
				}

				for (int& stamp: timestamps) {
					stamp = std::clamp(stamp, 0, static_cast<int>(t.size() - 1));
				}

				if (always_reset_axis) {
					solution_limits.assign(pde->start, pde->end, min_u, max_u, 0.05);
					error_limits.assign(0, t_end, 0, max_error, 0.05);

				}
			}

			ImGui::Separator();
			ImGui::Text("Timestamps");

			if (ImGui::Button("Add")) {
				timestamps.push_back(0);
			}
			ImGui::SameLine();
			if (ImGui::Button("Uniform")) {
				for (size_t i = 0; i < timestamps.size(); ++i) {
					timestamps[i] = static_cast<int>((t.size() - 1) * (double)i / (timestamps.size() - 1));
				}
			}

			if (!t.empty()) {
				for (size_t i = 0; i < timestamps.size(); ++i) {
					ImGui::BeginGroup();
					ImGui::SliderInt(format("##%d", i).get(), &timestamps[i], 0, t.size()-1);
					ImGui::SameLine();
					ImGui::Text(format("t = %.3f", t[timestamps[i]]).get());
					ImGui::SameLine();
					ImGui::PushID(i);
					if (ImGui::SmallButton("×")) {
						timestamps.erase(timestamps.begin() + i);
						i--;
					}
					ImGui::PopID();
					ImGui::EndGroup();
				}
			}

			ImGui::End();
		}

		ShowSolutionWindow("Explicit", x, t, u_expl, error_expl, timestamps, u_correct, show_correct_solution, solution_limits, error_limits, window_flags);
		ShowSolutionWindow("Implicit", x, t, u_impl, error_impl, timestamps, u_correct, show_correct_solution, solution_limits, error_limits, window_flags);
		ShowSolutionWindow("Combined", x, t, u_comb, error_comb, timestamps, u_correct, show_correct_solution, solution_limits, error_limits, window_flags);

		ImGui::Render();
		int display_w, display_h;
		glfwGetFramebufferSize(window, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
		glClear(GL_COLOR_BUFFER_BIT);
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImPlot::DestroyContext();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}
