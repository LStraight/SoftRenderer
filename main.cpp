#include "include/renderer.h"

int main() {
	Renderer window(800, 600);

	const int VC = 0;

	struct { vec4f pos; vec4f color; } vs_input[3] = {
		{ {  0.0,  0.7, 0.90, 1}, {1, 0, 0, 1} },
		{ { -0.6, -0.2, 0.01, 1}, {0, 1, 0, 1} },
		{ { +0.6, -0.2, 0.01, 1}, {0, 0, 1, 1} },
	};//图像中心取作原点（0，0）

	//index是三角形顶点序号，output是顶点上下文
	window.SetVertexShader([&](int index, ShaderContext& output)->vec4f {
		output.varying_vec4f[VC] = vs_input[index].color;
		return vs_input[index].pos;
		});

	window.SetPixelShader([&](ShaderContext& input)->vec4f {
		return input.varying_vec4f[VC];
		});
	window.DrawPrimitive();
	window.SaveFile("output.bmp");

	system("mspaint.exe output.bmp");

	return 0;
}