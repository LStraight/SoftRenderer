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
	//��������ͼԪ
	bool DrawPrimitive();

protected:
	struct Vertex {
		ShaderContext context;	//������
		float rhw;				//w�ĵ���
		vec4f pos;				//�ռ�����
		vec2f spf;				//������Ļ����
		vec2i spi;				//������Ļ����
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

	//�������㣬��ʼ��
	for (int k = 0; k < 3; ++k) {
		Vertex& vertex = m_vertex[k];

		vertex.context.varying_float.clear();
		vertex.context.varying_vec2f.clear();
		vertex.context.varying_vec3f.clear();
		vertex.context.varying_vec4f.clear();

		//��ȡ����λ����Ϣ
		vertex.pos = m_vertext_shader(k, vertex.context);
		float w = vertex.pos.w;

		//x,y,z��ȡֵ��Χ��[0,1)����w���߽��ж�������w��1
		//��һ���㳬���߽����ж�����ʧ�ܣ�����Ӧ�ý�Խ���߽�����������и�γ�����������������
		if (w == 0.f) return false;
		if (vertex.pos.x < -w || vertex.pos.x > w) return false;
		if (vertex.pos.y < -w || vertex.pos.y > w) return false;
		if (vertex.pos.z < 0.0f || vertex.pos.z > w) return false;

		//�����һ��
		vertex.rhw = 1.0f / w;
		vertex.pos *= vertex.rhw;
		vertex.spf.x = (vertex.pos.x + 1.0f) * fb_width * 0.5f;		//����������ת��Ϊ���ظ�������
		vertex.spf.y = (1.0f - vertex.pos.y) * fb_height * 0.5f;	
		vertex.spi.x = static_cast<int>(vertex.spf.x + 0.5f);		//ƫ��0.5�����������Ķ��룬ȡ��
		vertex.spi.y = static_cast<int>(vertex.spf.y + 0.5f);

		//����������
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

	//���������εıߣ��˴�Ӧ�п����
	if (render_frame) {
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[2].spi.x, m_vertex[2].spi.y);
		DrawLine(m_vertex[2].spi.x, m_vertex[2].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
	}

	//����������˳�
	if (render_pixel == false) return false;

	vec4f v01 = m_vertex[1].pos - m_vertex[0].pos;
	vec4f v02 = m_vertex[2].pos - m_vertex[0].pos;
	vec4f normal = vector_cross(v01, v02);

	//ʹ��vtx���ʶ��㣬��Ϊ��Ҫ��������ƽ��ķ����������˳��
	Vertex* vtx[3] = { &m_vertex[0], &m_vertex[1], &m_vertex[2] };

	//�涨�����η�������Ļ���
	if (normal.z > 0.0f) {
		vtx[1] = &m_vertex[2];
		vtx[2] = &m_vertex[1];
	}
	//���㹲�ߣ�����ʧ��
	else if (normal.z == 0.0f) {
		return false;
	}

	vec2i p0 = vtx[0]->spi;
	vec2i p1 = vtx[1]->spi;
	vec2i p2 = vtx[2]->spi;
	//�������
	float s = Abs(vector_cross(p1 - p0, p2 - p0));
	//���Ϊ0���˳�
	if (s <= 0) return false;

	//�ж������εıߵ�λ�ã�ֻ�����ϵı���Ҫ�����ɫ
	bool UpLeft01 = IsUpLeft(p0, p1);
	bool UpLeft12 = IsUpLeft(p1, p2);
	bool UpLeft20 = IsUpLeft(p2, p0);

	//������Ӿ����������ص㣬���������ڲ�����в�ֵ
	for (int cy = rect_min_y; cy <= rect_max_y; ++cy) {
		for (int cx = rect_min_x; cx <= rect_max_x; ++cx) {
			//����ƫ��0.5���ŵ���������
			vec2f px = { (float)cx + 0.5f, (float)cy + 0.5f };

			//��Ե����
			int E01 = -(cx - p0.x) * (p1.y - p0.y) + (cy - p0.y) * (p1.x - p0.x);
			int E12 = -(cx - p1.x) * (p2.y - p1.y) + (cy - p1.y) * (p2.x - p1.x);
			int E20 = -(cx - p2.x) * (p0.y - p2.y) + (cy - p2.y) * (p0.x - p2.x);

			//���ڱ�Ե�ϣ����ж��Ƿ�Ϊ���ϵıߣ����������䣬����continue
			if (E01 < (UpLeft01 ? 0 : 1)) continue;   
			if (E12 < (UpLeft12 ? 0 : 1)) continue;   
			if (E20 < (UpLeft20 ? 0 : 1)) continue;   

			//�˵㵽��ǰ�������
			vec2f s0 = vtx[0]->spf - px;
			vec2f s1 = vtx[1]->spf - px;
			vec2f s2 = vtx[2]->spf - px;

			float a = Abs(vector_cross(s1, s2));   
			float b = Abs(vector_cross(s2, s0));   
			float c = Abs(vector_cross(s0, s1));   
			float s = a + b + c;
			if (s == 0.0f) continue;

			//�ȱȷ�����ʹa+b+c = 1
			a = a * (1.0f / s);
			b = b * (1.0f / s);
			c = c * (1.0f / s);

			//�Ե�ǰ���rhw���������˵�����ֵ
			float rhw = vtx[0]->rhw * a + vtx[1]->rhw * b + vtx[2]->rhw * c;

			//��Ȳ���
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

			//Ϊ��ǰ��д��������context
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

	//���»������Ǳߣ����⸲��
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

	// �����ʼ��
	for (int k = 0; k < 3; k++) {
		Vertex& vertex = m_vertex[k];

		// ��������� varying �б�
		vertex.context.varying_float.clear();
		vertex.context.varying_vec2f.clear();
		vertex.context.varying_vec3f.clear();
		vertex.context.varying_vec4f.clear();

		// ���ж�����ɫ���򣬷��ض�������
		vertex.pos = m_vertext_shader(k, vertex.context);

		// �򵥲ü����κ�һ�����㳬�� CVV ���޳�
		float w = vertex.pos.w;

		// ����ͼ�򵥣���һ����Խ�磬����������������Σ�����ϸ��������
		// ���Խ���˾�����οռ��ڽ��вü������Ϊ 0-2 ��������Ȼ�����
		if (w == 0.0f) return false;
		if (vertex.pos.z < 0.0f || vertex.pos.z > w) return false;
		if (vertex.pos.x < -w || vertex.pos.x > w) return false;
		if (vertex.pos.y < -w || vertex.pos.y > w) return false;

		// ���� w �ĵ�����Reciprocal of the Homogeneous W 
		vertex.rhw = 1.0f / w;

		// �������ռ� /w ��һ������λ��� cvv
		vertex.pos *= vertex.rhw;

		// ������Ļ����
		vertex.spf.x = (vertex.pos.x + 1.0f) * fb_width * 0.5f;
		vertex.spf.y = (1.0f - vertex.pos.y) * fb_height * 0.5f;

		// ������Ļ���꣺�� 0.5 ��ƫ��ȡ��Ļ���ط������Ķ���
		vertex.spi.x = (int)(vertex.spf.x + 0.5f);
		vertex.spi.y = (int)(vertex.spf.y + 0.5f);

		// ������Ӿ��η�Χ
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

	// �����߿�
	if (render_frame) {
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[2].spi.x, m_vertex[2].spi.y);
		DrawLine(m_vertex[2].spi.x, m_vertex[2].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
	}

	// �����������ؾ��˳�
	if (render_pixel == false) return false;

	// �ж������γ���
	vec4f v01 = m_vertex[1].pos - m_vertex[0].pos;
	vec4f v02 = m_vertex[2].pos - m_vertex[0].pos;
	vec4f normal = vector_cross(v01, v02);

	// ʹ�� vtx �����������㣬����ֱ���� _vertex ���ʣ���Ϊ���ܻ����˳��
	Vertex* vtx[3] = { &m_vertex[0], &m_vertex[1], &m_vertex[2] };

	// ��������ӵ㣬�򽻻����㣬��֤ edge equation �жϵķ���Ϊ��
	if (normal.z > 0.0f) {
		vtx[1] = &m_vertex[2];
		vtx[2] = &m_vertex[1];
	}
	else if (normal.z == 0.0f) {
		return false;
	}

	// ���������˵�λ��
	vec2i p0 = vtx[0]->spi;
	vec2i p1 = vtx[1]->spi;
	vec2i p2 = vtx[2]->spi;

	// ���������Ϊ����˳�
	float s = Abs(vector_cross(p1 - p0, p2 - p0));
	if (s <= 0) return false;

	// ���������ʱ�����������ı��ϵĵ���Ҫ�������ҷ����·����ϵĵ㲻����
	// ���ж��Ƿ��� TopLeft���жϳ��������·� Edge Equation һ�����
	bool TopLeft01 = IsUpLeft(p0, p1);
	bool TopLeft12 = IsUpLeft(p1, p2);
	bool TopLeft20 = IsUpLeft(p2, p0);

	// ������������Ӿ��ε����е�
	for (int cy = rect_min_y; cy <= rect_max_y; cy++) {
		for (int cx = rect_min_x; cx <= rect_max_x; cx++) {
			vec2f px = { (float)cx + 0.5f, (float)cy + 0.5f };

			// Edge Equation
			// ʹ���������⸡����ͬʱ��Ϊ������ϵ�����Է���ȡ��
			int E01 = -(cx - p0.x) * (p1.y - p0.y) + (cy - p0.y) * (p1.x - p0.x);
			int E12 = -(cx - p1.x) * (p2.y - p1.y) + (cy - p1.y) * (p2.x - p1.x);
			int E20 = -(cx - p2.x) * (p0.y - p2.y) + (cy - p2.y) * (p0.x - p2.x);


			// ��������ϱߣ��� E >= 0 �жϺϷ���������±߾��� E > 0 �жϺϷ�
			// ����ͨ������һ����� 1 ������ < 0 �� <= 0 ��һ��ʽ�ӱ��
			if (E01 < (TopLeft01 ? 0 : 1)) continue;   // �ڵ�һ���ߺ���
			if (E12 < (TopLeft12 ? 0 : 1)) continue;   // �ڵڶ����ߺ���
			if (E20 < (TopLeft20 ? 0 : 1)) continue;   // �ڵ������ߺ���

			// �����˵㵽��ǰ���ʸ��
			vec2f s0 = vtx[0]->spf - px;
			vec2f s1 = vtx[1]->spf - px;
			vec2f s2 = vtx[2]->spf - px;

			// ��������ϵ�������ڲ������������ a / b / c
			float a = Abs(vector_cross(s1, s2));    // �������� Px-P1-P2 ���
			float b = Abs(vector_cross(s2, s0));    // �������� Px-P2-P0 ���
			float c = Abs(vector_cross(s0, s1));    // �������� Px-P0-P1 ���
			float s = a + b + c;                    // �������� P0-P1-P2 ���

			if (s == 0.0f) continue;

			// ������������Ա�֤��a + b + c = 1������������ֵϵ��
			a = a * (1.0f / s);
			b = b * (1.0f / s);
			c = c * (1.0f / s);

			// ���㵱ǰ��� 1/w���� 1/w ����Ļ�ռ�����Թ�ϵ����ֱ�����Ĳ�ֵ
			float rhw = vtx[0]->rhw * a + vtx[1]->rhw * b + vtx[2]->rhw * c;

			// ������Ȳ���
			if (rhw <  depth_buffer[cy][cx]) continue;
			 depth_buffer[cy][cx] = rhw;   // ��¼ 1/w ����Ȼ���

			// ��ԭ��ǰ���ص� w
			float w = 1.0f / ((rhw != 0.0f) ? rhw : 1.0f);

			// �������������ֵ varying ��ϵ��
			// �ȳ��Ը��Զ���� w Ȼ�������Ļ�ռ��ֵȻ���ٳ��Ե�ǰ w
			float c0 = vtx[0]->rhw * a * w;
			float c1 = vtx[1]->rhw * b * w;
			float c2 = vtx[2]->rhw * c * w;

			// ׼��Ϊ��ǰ���صĸ��� varying ���в�ֵ
			ShaderContext input;

			ShaderContext& i0 = vtx[0]->context;
			ShaderContext& i1 = vtx[1]->context;
			ShaderContext& i2 = vtx[2]->context;

			// ��ֵ���� varying
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

			// ִ��������ɫ��
			vec4f color = { 0.0f, 0.0f, 0.0f, 0.0f };

			if (m_pixel_shader != NULL) {
				color = m_pixel_shader(input);
			}

			// ���Ƶ� framebuffer �ϣ�������Լ��жϣ���� PS ���ص���ɫ alpha ����
			// С�ڵ�������������ƣ����������Ļ�Ҫ��ǰ��ĸ�����Ȼ���Ĵ���Ų������
			// ֻ����Ҫ��Ⱦ��ʱ��Ÿ�����ȡ�
			frame_buffer->SetPixel(cx, cy, color);
		}
	}

	// �����߿��ٻ�һ�α��⸲��
	if (render_frame) {
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
		DrawLine(m_vertex[0].spi.x, m_vertex[0].spi.y, m_vertex[2].spi.x, m_vertex[2].spi.y);
		DrawLine(m_vertex[2].spi.x, m_vertex[2].spi.y, m_vertex[1].spi.x, m_vertex[1].spi.y);
	}

	return true;
}

*/