#pragma once

#include <map>
#include <functional>

#include "vec.h"
#include "matrix.h"
#include "mathmatics.h"
#include "Bitmap.h"

struct ShaderContext {
	std::map<int, float> varying_float;
	std::map<int, vec2f> varying_vec2f;
	std::map<int, vec3f> varying_vec3f;
	std::map<int, vec4f> varying_vec4f;
};

typedef std::function<vec4f(int index, ShaderContext& output)> VertexShader;
typedef std::function<vec4f(ShaderContext& input)> PixelShader;

class Renderer {
public:
	Renderer();
	Renderer(int width, int height);
	virtual ~Renderer() { Reset(); }

public:
	void Reset();
	void Init(int width, int height);
	void Clear();
	void SetVertexShader(VertexShader vs);
	void SetPixelShader(PixelShader ps);
	void SaveFile(const char* filename);
	void SetBackgroundColor(uint32_t color);
	void SetForegroundColor(uint32_t color);
	void SetPixel(int x, int y, uint32_t cc);
	void SetPixel(int x, int y, const vec4f cc);
	void SetPixel(int x, int y, const vec3f cc);
	void DrawLine(int x1, int y1, int x2, int y2);
	void SetRenderState(bool frame, bool pixel);
	bool IsUpLeft(const vec2i& a, const vec2i& b);
	//绘制三角图元
	bool DrawPrimitive();

protected:
	struct Vertex {
		ShaderContext context;	//上下文
		float rhw;				//w的倒数
		vec4f pos;				//空间坐标
		vec2f spf;				//浮点屏幕坐标
		vec2i spi;				//整数屏幕坐标
	};

	Bitmap* frame_buffer;
	float** depth_buffer;
	int fb_width;
	int fb_height;
	uint32_t foreground_color;
	uint32_t background_color;

	Vertex m_vertex[3];
	int rect_min_x;
	int rect_min_y;
	int rect_max_x;
	int rect_max_y;

	bool render_frame;
	bool render_pixel;

	VertexShader m_vertext_shader;
	PixelShader m_pixel_shader;
};

Renderer::Renderer() {
	frame_buffer = nullptr;
	depth_buffer = nullptr;
	render_frame = false;
	render_pixel = true;
}

Renderer::Renderer(int width, int height) {
	frame_buffer = nullptr;
	depth_buffer = nullptr;
	render_frame = false;
	render_pixel = true;
	Init(width, height);
}

void Renderer::Reset() {
	m_vertext_shader = nullptr;
	m_pixel_shader = nullptr;

	if (frame_buffer)
		delete frame_buffer;
	frame_buffer = nullptr;

	if (depth_buffer) {
		for (int i = 0; i < fb_height; ++i) {
			if (depth_buffer[i])
				delete[] depth_buffer[i];
			depth_buffer[i] = nullptr;
		}
		delete[] depth_buffer;
		depth_buffer = nullptr;
	}
	foreground_color = 0xffffffff;
	background_color = 0xff191970;
}
void Renderer::Init(int width, int height) {
	Reset();
	frame_buffer = new Bitmap(width, height);
	fb_width = width;
	fb_height = height;
	depth_buffer = new float* [height];
	for (int i = 0; i < height; ++i) {
		depth_buffer[i] = new float[width];
	}
	Clear();
}

void Renderer::Clear() {
	if (frame_buffer)
		frame_buffer->Fill(background_color);
	if(depth_buffer)
		for (int i = 0; i < fb_height; ++i) {
			for (int j = 0; j < fb_width; ++j) {
				depth_buffer[i][j] = 0.0f;
			}
		}
}

void Renderer::SetVertexShader(VertexShader vs) {
	m_vertext_shader = vs;
}

void Renderer::SetPixelShader(PixelShader ps) {
	m_pixel_shader = ps;
}

void Renderer::SaveFile(const char* filename) {
	if (frame_buffer)
		frame_buffer->SaveFile(filename);
}

void Renderer::SetBackgroundColor(uint32_t color) {
	background_color = color;
}

void Renderer::SetForegroundColor(uint32_t color) {
	foreground_color = color;
}

void Renderer::SetPixel(int x, int y, uint32_t cc) {
	if (frame_buffer)
		frame_buffer->SetPixel(x, y, cc);
}

void Renderer::SetPixel(int x, int y, const vec4f cc) {
	SetPixel(x, y, vector_to_color(cc));
}

void Renderer::SetPixel(int x, int y, const vec3f cc) {
	SetPixel(x, y, vector_to_color(cc));
}

void Renderer::DrawLine(int x1, int y1, int x2, int y2) {
	if (frame_buffer)
		frame_buffer->DrawLine(x1, y1, x2, y2, foreground_color);
}

void Renderer::SetRenderState(bool frame, bool pixel) {
	render_frame = frame;
	render_pixel = pixel;
}

bool Renderer::IsUpLeft(const vec2i& a, const vec2i& b) {
	return ((a.y == b.y) && (a.x < b.x)) || (a.y > b.y);
}

bool Renderer::DrawPrimitive() {
	if (!frame_buffer || !m_vertext_shader)
		std::cerr << "draw failed" << std::endl;
		return false;

	//遍历顶点，初始化
	for (int k = 0; k < 3; ++k) {
		Vertex& vertex = m_vertex[k];

		vertex.context.varying_float.clear();
		vertex.context.varying_vec2f.clear();
		vertex.context.varying_vec3f.clear();
		vertex.context.varying_vec4f.clear();

		//读取顶点位置信息
		vertex.pos = m_vertext_shader(k, vertex.context);
		float w = vertex.pos.w;

		//x,y,z的取值范围在[0,1)，用w做边界判定，这里w是1
		//有一个点超出边界这判定绘制失败，正常应该将越出边界的三角形做切割，形成两个完整的三角形
		if (w == 0.f) return false;
		if (vertex.pos.x < -w || vertex.pos.x > w) return false;
		if (vertex.pos.y < -w || vertex.pos.y > w) return false;
		if (vertex.pos.z < 0.0f || vertex.pos.z > w) return false;

		//坐标归一化
		vertex.rhw = 1.0f / w;
		vertex.pos *= vertex.rhw;
		vertex.spf.x = (vertex.pos.x + 1.0f) * fb_width * 0.5f;		//将比例坐标转化为像素浮点坐标
		vertex.spf.y = (1.0f - vertex.pos.y) * fb_height * 0.5f;	
		vertex.spi.x = static_cast<int>(vertex.spf.x + 0.5f);		//偏置0.5，与像素中心对齐，取整
		vertex.spi.y = static_cast<int>(vertex.spf.y + 0.5f);

		//外界矩形坐标
		if (k == 0) {
			rect_min_x = rect_max_x = Between(0, fb_width - 1, vertex.spi.x);
			rect_min_y = rect_max_y = Between(0, fb_height - 1, vertex.spi.y);
		}
		else {
			rect_min_x = Between(0, fb_width - 1, Min(rect_min_x, vertex.spi.x));
			rect_max_x = Between(0, fb_width - 1, Max(rect_max_x, vertex.spi.x));
			rect_min_y = Between(0, fb_height- 1, Min(rect_min_y, vertex.spi.y));
			rect_min_x = Between(0, fb_height- 1, Max(rect_max_y, vertex.spi.y));
		}
	}

	//绘制三角形的边（此处应有抗锯齿
	if (render_frame) {
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[2].spi.x, m_vertex[2].spi.y);
		DrawLine(m_vertex[2].spi.x, m_vertex[2].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
	}

	//无需填充则退出
	if (render_pixel == false) return false;

	vec4f v01 = m_vertex[1].pos - m_vertex[0].pos;
	vec4f v02 = m_vertex[2].pos - m_vertex[0].pos;
	vec4f normal = vector_cross(v01, v02);

	//使用vtx访问顶点，因为需要根据三角平面的法向调整顶点顺序
	Vertex* vtx[3] = { &m_vertex[0], &m_vertex[1], &m_vertex[2] };

	//规定三角形法向朝向屏幕里侧
	if (normal.z > 0.0f) {
		vtx[1] = &m_vertex[2];
		vtx[2] = &m_vertex[1];
	}
	//三点共线，绘制失败
	else if (normal.z == 0.0f) {
		return false;
	}

	vec2i p0 = vtx[0]->spi;
	vec2i p1 = vtx[1]->spi;
	vec2i p2 = vtx[2]->spi;
	//计算面积
	float s = Abs(vector_cross(p1 - p0, p2 - p0));
	//面积为0，退出
	if (s <= 0) return false;

	//判断三角形的边的位置，只有左、上的边需要填充颜色
	bool UpLeft01 = IsUpLeft(p0, p1);
	bool UpLeft12 = IsUpLeft(p1, p2);
	bool UpLeft20 = IsUpLeft(p2, p0);

	//遍历外接矩阵所有像素点，对三角形内部点进行插值
	for (int cy = rect_min_y; cy <= rect_max_y; ++cy) {
		for (int cx = rect_min_x; cx <= rect_max_x; ++cx) {
			//坐标偏置0.5，放到像素中心
			vec2f px = { (float)cx + 0.5f, (float)cy + 0.5f };

			//边缘方程
			int E01 = -(cx - p0.x) * (p1.y - p0.y) + (cy - p0.y) * (p1.x - p0.x);
			int E12 = -(cx - p1.x) * (p2.y - p1.y) + (cy - p1.y) * (p2.x - p1.x);
			int E20 = -(cx - p2.x) * (p0.y - p2.y) + (cy - p2.y) * (p0.x - p2.x);

			//点在边缘上，则判断是否为左、上的边，是则继续填充，否则continue
			if (E01 < (UpLeft01 ? 0 : 1)) continue;   
			if (E12 < (UpLeft12 ? 0 : 1)) continue;   
			if (E20 < (UpLeft20 ? 0 : 1)) continue;   

			//端点到当前点的向量
			vec2f s0 = vtx[0]->spf - px;
			vec2f s1 = vtx[1]->spf - px;
			vec2f s2 = vtx[2]->spf - px;

			float a = Abs(vector_cross(s1, s2));   
			float b = Abs(vector_cross(s2, s0));   
			float c = Abs(vector_cross(s0, s1));   
			float s = a + b + c;
			if (s == 0.0f) continue;

			//等比放缩，使a+b+c = 1
			a = a * (1.0f / s);
			b = b * (1.0f / s);
			c = c * (1.0f / s);

			//对当前点的rhw关于三个端点做插值
			float rhw = vtx[0]->rhw * a + vtx[1]->rhw * b + vtx[2]->rhw * c;

			//深度测试
			if (rhw < depth_buffer[cy][cx]) continue;
			depth_buffer[cy][cx] = rhw;

			float w = 1.0f / ((rhw != 0.0f) ? rhw : 1.0f);
			float c0 = vtx[0]->rhw * a * w;
			float c1 = vtx[1]->rhw * b * w;
			float c2 = vtx[2]->rhw * c * w;

			ShaderContext input;
			ShaderContext& i0 = vtx[0]->context;
			ShaderContext& i1 = vtx[1]->context;
			ShaderContext& i2 = vtx[2]->context;

			//为当前点写入上下文context
			for (auto const& it : i0.varying_float) {
				int key = it.first;
				float f0 = i0.varying_float[key];
				float f1 = i1.varying_float[key];
				float f2 = i2.varying_float[key];
				input.varying_float[key] = c0 * f0 + c1 * f1 + c2 * f2;
			}
			for (auto const& it : i0.varying_vec2f) {
				int key = it.first;
				const vec2f& f0 = i0.varying_vec2f[key];
				const vec2f& f1 = i1.varying_vec2f[key];
				const vec2f& f2 = i2.varying_vec2f[key];
				input.varying_vec2f[key] = c0 * f0 + c1 * f1 + c2 * f2;
			}
			for (auto const& it : i0.varying_vec3f) {
				int key = it.first;
				const vec3f& f0 = i0.varying_vec3f[key];
				const vec3f& f1 = i1.varying_vec3f[key];
				const vec3f& f2 = i2.varying_vec3f[key];
				input.varying_vec3f[key] = c0 * f0 + c1 * f1 + c2 * f2;
			}
			for (auto const& it : i0.varying_vec4f) {
				int key = it.first;
				const vec4f& f0 = i0.varying_vec4f[key];
				const vec4f& f1 = i1.varying_vec4f[key];
				const vec4f& f2 = i2.varying_vec4f[key];
				input.varying_vec4f[key] = c0 * f0 + c1 * f1 + c2 * f2;
			}

			vec4f color = { 0.0f, 0.0f, 0.0f, 0.0f };
			if (m_pixel_shader != NULL) {
				color = m_pixel_shader(input);
			}
			frame_buffer->SetPixel(cx, cy, color);
		}
	}

	//重新绘制三角边，避免覆盖
	if (render_frame) {
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[2].spi.x, m_vertex[2].spi.y);
		DrawLine(m_vertex[2].spi.x, m_vertex[2].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
	}
	return true;
}



/*
bool Renderer::DrawPrimitive() {
	if (frame_buffer == NULL || m_vertext_shader == NULL)
		return false;

	// 顶点初始化
	for (int k = 0; k < 3; k++) {
		Vertex& vertex = m_vertex[k];

		// 清空上下文 varying 列表
		vertex.context.varying_float.clear();
		vertex.context.varying_vec2f.clear();
		vertex.context.varying_vec3f.clear();
		vertex.context.varying_vec4f.clear();

		// 运行顶点着色程序，返回顶点坐标
		vertex.pos = m_vertext_shader(k, vertex.context);

		// 简单裁剪，任何一个顶点超过 CVV 就剔除
		float w = vertex.pos.w;

		// 这里图简单，当一个点越界，立马放弃整个三角形，更精细的做法是
		// 如果越界了就在齐次空间内进行裁剪，拆分为 0-2 个三角形然后继续
		if (w == 0.0f) return false;
		if (vertex.pos.z < 0.0f || vertex.pos.z > w) return false;
		if (vertex.pos.x < -w || vertex.pos.x > w) return false;
		if (vertex.pos.y < -w || vertex.pos.y > w) return false;

		// 计算 w 的倒数：Reciprocal of the Homogeneous W 
		vertex.rhw = 1.0f / w;

		// 齐次坐标空间 /w 归一化到单位体积 cvv
		vertex.pos *= vertex.rhw;

		// 计算屏幕坐标
		vertex.spf.x = (vertex.pos.x + 1.0f) * fb_width * 0.5f;
		vertex.spf.y = (1.0f - vertex.pos.y) * fb_height * 0.5f;

		// 整数屏幕坐标：加 0.5 的偏移取屏幕像素方格中心对齐
		vertex.spi.x = (int)(vertex.spf.x + 0.5f);
		vertex.spi.y = (int)(vertex.spf.y + 0.5f);

		// 更新外接矩形范围
		if (k == 0) {
			rect_min_x = rect_max_x = Between(0, fb_width - 1, vertex.spi.x);
			rect_min_y = rect_max_y = Between(0, fb_height - 1, vertex.spi.y);
		}
		else {
			rect_min_x = Between(0, fb_width - 1,  Min(rect_min_x, vertex.spi.x));
			rect_max_x = Between(0, fb_width - 1,  Max(rect_max_x, vertex.spi.x));
			rect_min_y = Between(0, fb_height - 1, Min(rect_min_y, vertex.spi.y));
			rect_max_y = Between(0, fb_height - 1, Max(rect_max_y, vertex.spi.y));
		}
	}

	// 绘制线框
	if (render_frame) {
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[2].spi.x, m_vertex[2].spi.y);
		DrawLine(m_vertex[2].spi.x, m_vertex[2].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
	}

	// 如果不填充像素就退出
	if (render_pixel == false) return false;

	// 判断三角形朝向
	vec4f v01 = m_vertex[1].pos - m_vertex[0].pos;
	vec4f v02 = m_vertex[2].pos - m_vertex[0].pos;
	vec4f normal = vector_cross(v01, v02);

	// 使用 vtx 访问三个顶点，而不直接用 _vertex 访问，因为可能会调整顺序
	Vertex* vtx[3] = { &m_vertex[0], &m_vertex[1], &m_vertex[2] };

	// 如果背向视点，则交换顶点，保证 edge equation 判断的符号为正
	if (normal.z > 0.0f) {
		vtx[1] = &m_vertex[2];
		vtx[2] = &m_vertex[1];
	}
	else if (normal.z == 0.0f) {
		return false;
	}

	// 保存三个端点位置
	vec2i p0 = vtx[0]->spi;
	vec2i p1 = vtx[1]->spi;
	vec2i p2 = vtx[2]->spi;

	// 计算面积，为零就退出
	float s = Abs(vector_cross(p1 - p0, p2 - p0));
	if (s <= 0) return false;

	// 三角形填充时，左面和上面的边上的点需要包括，右方和下方边上的点不包括
	// 先判断是否是 TopLeft，判断出来后会和下方 Edge Equation 一起决策
	bool TopLeft01 = IsUpLeft(p0, p1);
	bool TopLeft12 = IsUpLeft(p1, p2);
	bool TopLeft20 = IsUpLeft(p2, p0);

	// 迭代三角形外接矩形的所有点
	for (int cy = rect_min_y; cy <= rect_max_y; cy++) {
		for (int cx = rect_min_x; cx <= rect_max_x; cx++) {
			vec2f px = { (float)cx + 0.5f, (float)cy + 0.5f };

			// Edge Equation
			// 使用整数避免浮点误差，同时因为是左手系，所以符号取反
			int E01 = -(cx - p0.x) * (p1.y - p0.y) + (cy - p0.y) * (p1.x - p0.x);
			int E12 = -(cx - p1.x) * (p2.y - p1.y) + (cy - p1.y) * (p2.x - p1.x);
			int E20 = -(cx - p2.x) * (p0.y - p2.y) + (cy - p2.y) * (p0.x - p2.x);


			// 如果是左上边，用 E >= 0 判断合法，如果右下边就用 E > 0 判断合法
			// 这里通过引入一个误差 1 ，来将 < 0 和 <= 0 用一个式子表达
			if (E01 < (TopLeft01 ? 0 : 1)) continue;   // 在第一条边后面
			if (E12 < (TopLeft12 ? 0 : 1)) continue;   // 在第二条边后面
			if (E20 < (TopLeft20 ? 0 : 1)) continue;   // 在第三条边后面

			// 三个端点到当前点的矢量
			vec2f s0 = vtx[0]->spf - px;
			vec2f s1 = vtx[1]->spf - px;
			vec2f s2 = vtx[2]->spf - px;

			// 重心坐标系：计算内部子三角形面积 a / b / c
			float a = Abs(vector_cross(s1, s2));    // 子三角形 Px-P1-P2 面积
			float b = Abs(vector_cross(s2, s0));    // 子三角形 Px-P2-P0 面积
			float c = Abs(vector_cross(s0, s1));    // 子三角形 Px-P0-P1 面积
			float s = a + b + c;                    // 大三角形 P0-P1-P2 面积

			if (s == 0.0f) continue;

			// 除以总面积，以保证：a + b + c = 1，方便用作插值系数
			a = a * (1.0f / s);
			b = b * (1.0f / s);
			c = c * (1.0f / s);

			// 计算当前点的 1/w，因 1/w 和屏幕空间呈线性关系，故直接重心插值
			float rhw = vtx[0]->rhw * a + vtx[1]->rhw * b + vtx[2]->rhw * c;

			// 进行深度测试
			if (rhw <  depth_buffer[cy][cx]) continue;
			 depth_buffer[cy][cx] = rhw;   // 记录 1/w 到深度缓存

			// 还原当前像素的 w
			float w = 1.0f / ((rhw != 0.0f) ? rhw : 1.0f);

			// 计算三个顶点插值 varying 的系数
			// 先除以各自顶点的 w 然后进行屏幕空间插值然后再乘以当前 w
			float c0 = vtx[0]->rhw * a * w;
			float c1 = vtx[1]->rhw * b * w;
			float c2 = vtx[2]->rhw * c * w;

			// 准备为当前像素的各项 varying 进行插值
			ShaderContext input;

			ShaderContext& i0 = vtx[0]->context;
			ShaderContext& i1 = vtx[1]->context;
			ShaderContext& i2 = vtx[2]->context;

			// 插值各项 varying
			for (auto const& it : i0.varying_float) {
				int key = it.first;
				float f0 = i0.varying_float[key];
				float f1 = i1.varying_float[key];
				float f2 = i2.varying_float[key];
				input.varying_float[key] = c0 * f0 + c1 * f1 + c2 * f2;
			}

			for (auto const& it : i0.varying_vec2f) {
				int key = it.first;
				const vec2f& f0 = i0.varying_vec2f[key];
				const vec2f& f1 = i1.varying_vec2f[key];
				const vec2f& f2 = i2.varying_vec2f[key];
				input.varying_vec2f[key] = c0 * f0 + c1 * f1 + c2 * f2;
			}

			for (auto const& it : i0.varying_vec3f) {
				int key = it.first;
				const vec3f& f0 = i0.varying_vec3f[key];
				const vec3f& f1 = i1.varying_vec3f[key];
				const vec3f& f2 = i2.varying_vec3f[key];
				input.varying_vec3f[key] = c0 * f0 + c1 * f1 + c2 * f2;
			}

			for (auto const& it : i0.varying_vec4f) {
				int key = it.first;
				const vec4f& f0 = i0.varying_vec4f[key];
				const vec4f& f1 = i1.varying_vec4f[key];
				const vec4f& f2 = i2.varying_vec4f[key];
				input.varying_vec4f[key] = c0 * f0 + c1 * f1 + c2 * f2;
			}

			// 执行像素着色器
			vec4f color = { 0.0f, 0.0f, 0.0f, 0.0f };

			if (m_pixel_shader != NULL) {
				color = m_pixel_shader(input);
			}

			// 绘制到 framebuffer 上，这里可以加判断，如果 PS 返回的颜色 alpha 分量
			// 小于等于零则放弃绘制，不过这样的话要把前面的更新深度缓存的代码挪下来，
			// 只有需要渲染的时候才更新深度。
			frame_buffer->SetPixel(cx, cy, color);
		}
	}

	// 绘制线框，再画一次避免覆盖
	if (render_frame) {
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[2].spi.x, m_vertex[2].spi.y);
		DrawLine(m_vertex[2].spi.x, m_vertex[2].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
	}

	return true;
}

*/