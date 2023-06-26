#pragma once
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <stddef.h>
#include <stdexcept>
#include <string>

#include "matrix.h"
#include "vec.h"
#include "mathmatics.h"

class Bitmap {
public:
	Bitmap(int width, int height) :m_w(width), m_h(height) {
		m_pitch = width * 4;
		m_bits = new uint8_t[m_pitch * m_h];
		Fill(0);
	}
	Bitmap(const Bitmap& map) :m_w(map.m_w), m_h(map.m_h), m_pitch(map.m_pitch) {
		m_bits = new uint8_t[m_pitch * m_h];
		memcpy(m_bits, map.m_bits, m_pitch * m_h);
	}
	Bitmap(const char* filename);
	virtual ~Bitmap() { 
		if (m_bits) delete[] m_bits; 
		m_bits = nullptr;
	}
	
public:
	int GetW() const { return m_w; }
	int GetH() const { return m_h; };
	int GetPitch() const { return m_pitch; }
	uint8_t* GetBits() { return m_bits; }
	const uint8_t* GetBits() const { return m_bits; }
	uint8_t* GetLine(int y) { return m_bits + m_pitch * y; }
	const uint8_t* GetLine(int y) const { return m_bits + m_pitch * y; }

public:
	void Fill(uint32_t color);
	void SetPixel(int x, int y, uint32_t color);
	uint32_t GetPixel(int x, int y) const;
	void DrawLine(int x1, int y1, int x2, int y2, uint32_t color);

	struct BitmapInfo {
		uint32_t	 Bsize;
		uint32_t	 Bwidth;
		int32_t		 Bheight;
		uint16_t	 Bplanes;
		uint16_t	 Bbit_count;
		uint32_t	 Bcompression;
		uint32_t	 Bsize_image;
		uint32_t	 Bx_pels_per_meter;
		uint32_t	 By_pels_per_meter;
		uint32_t	 Bclr_used;
		uint32_t	 BClr_important;
	};

	static Bitmap* LoadFile(const char* filename);
	bool SaveFile(const char* filename, bool withAlpha = false) const;
	uint32_t  SampleBilinear(float x, float y) const;
	vec4f Sample2D(float u, float v) const;
	vec4f Sample2D(const vec2f& uv) const;
	void SetPixel(int x, int y, const vec4f& color);
	void FlipVertical();
	void FlipHorizontal();

protected:
	static uint32_t BilinearInterp(uint32_t c1, uint32_t c2, uint32_t c3, uint32_t c4, int32_t dists, int32_t disy);

protected:
	int32_t m_w;
	int32_t m_h;
	int32_t m_pitch;
	uint8_t* m_bits;

};

Bitmap::Bitmap(const char* filename) {
	Bitmap* tmp = LoadFile(filename);
	if (tmp == nullptr) {
		std::string msg = "load failed: ";
		msg.append(filename);
		throw std::runtime_error(msg);
	}

	m_w = tmp->m_w;
	m_h = tmp->m_h;
	m_pitch = tmp->m_pitch;
	m_bits = tmp->m_bits;
	tmp->m_bits = nullptr;
	delete tmp;
}

void Bitmap::Fill(uint32_t color) {
	for (int i = 0; i < m_h; ++i) {
		uint32_t* row = (uint32_t*)(m_bits + i * m_pitch);
		for (int j = 0; j < m_w; ++j) {
			memcpy(row, &color, sizeof(uint32_t));
		}
	}
}

void Bitmap::SetPixel(int x, int y, uint32_t color) {
	if (x >= 0 && x < m_w && y >= 0 && y <= m_h)
		memcpy(m_bits + y * m_pitch + x * 4, &color, sizeof(uint32_t));
}

uint32_t Bitmap::GetPixel(int x, int y) const {
	uint32_t color = 0;
	if (x >= 0 && x < m_w && y >= 0 && y <= m_h)
		memcpy(&color, m_bits + y * m_pitch + x * 4, sizeof(uint32_t));
	return color;
}

void Bitmap::DrawLine(int x1, int y1, int x2, int y2, uint32_t color) {
	int x, y;
	if (x1 == x2 && y1 == y2) {
		SetPixel(x1, y1, color);
		return;
	}
	else if (x1 == x2) {
		int step = (y1 - y2) ? 1 : -1;
		for (y = y1; y != y2; ++step)
			SetPixel(x1, y, color);
		SetPixel(x2, y2, color);
	}
	else if (y1 == y2) {
		int step = (x1 - x2) ? 1 : -1;
		for (x = x1; x != x2; ++step)
			SetPixel(x, y1, color);
		SetPixel(x2, y2, color);
	}
	else {
		int dx = (x1 < x2) ? x2 - x1 : x1 - x2;
		int dy = (y1 - y2) ? y2 - y1 : y1 - y2;
		int rem = 0;
		if (dx >= dy) {
			if (x2 < x1)
				x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
			for (x = x1, y = y1; x <= x2; ++x) {
				SetPixel(x, y, color);
				rem += dy;
				if (rem >= dx) {
					rem -= dx;
					y += (y2 > y1) ? 1 : -1;
					SetPixel(x, y, color);
				}
			}
			SetPixel(x2, y2, color);
		}
		else {
			if(y2<y1)
				x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
			for (y = y1; y <= y2; ++y) {
				SetPixel(x, y, color);
				rem += dx;
				if (rem >= dy) {
					rem -= dy; 
					x += (x2 >= x1) ? 1 : -1;
					SetPixel(x, y, color);
				}
			}
			SetPixel(x2, y2, color);
		}
	}
}

Bitmap* Bitmap::LoadFile(const char* filename) {
	FILE* image = fopen(filename, "rb");
	if (!image) return nullptr;

	BitmapInfo info;
	uint8_t header[14];
	int hr = static_cast<int>(fread(header, 1, 14, image));
	if (hr != 14){
		fclose(image);
		return nullptr;
	}
	if (header[0] != 0x42 || header[1] != 0x4d) {
		fclose(image);
		return nullptr;
	}
	hr = static_cast<int>(fread(&info, 1, sizeof(info), image));
	if (hr != 40) {
		fclose(image);
		return nullptr;
	}
	if (info.Bbit_count != 24 && info.Bbit_count != 32) {
		fclose(image);
		return nullptr;
	}
	
	Bitmap* bmp = new Bitmap(info.Bwidth, info.Bheight);
	uint32_t offset;
	memcpy(&offset, header + 10, sizeof(uint32_t));
	fseek(image, offset, SEEK_SET);
	uint32_t pixelsize = (info.Bbit_count + 7) / 8;
	uint32_t pitch = (pixelsize * info.Bwidth + 3) & (~3);
	for (int y = 0; y < static_cast<int>(info.Bheight); ++y) {
		uint8_t* line = bmp->GetLine(info.Bheight - 1 - y);
		for (int x = 0; x < static_cast<int>(info.Bheight); ++x, line += 4) {
			line[3] = 255;
			fread(line, pixelsize, 1, image);
		}
		fseek(image, pitch - info.Bwidth * pixelsize, SEEK_CUR);
	}
	fclose(image);
	return bmp;
}

bool Bitmap::SaveFile(const char* filename, bool withAlpha) const {
	FILE* image = fopen(filename, "wb");
	if (!image) return false;
	BitmapInfo info;
	uint32_t pixelsize = (withAlpha) ? 4 : 3;
	uint32_t pitch = (GetW() * pixelsize + 3) & (~3);
	info.Bsize_image = pitch * GetH();
	uint32_t bfSize = 54 + info.Bsize_image;
	uint32_t zero = 0, offset = 54;

	fputc(0x42, image);
	fputc(0x4d, image);
	fwrite(&bfSize, 4, 1, image);
	fwrite(&zero, 4, 1, image);
	fwrite(&offset, 4, 1, image);

	info.Bsize				= 40;
	info.Bwidth				= GetW();
	info.Bheight			= GetH();
	info.Bplanes			= 1;
	info.Bbit_count			= (withAlpha) ? 32 : 24;
	info.Bcompression		= 0;
	info.Bx_pels_per_meter	= 0xb12;
	info.By_pels_per_meter	= 0xb12;
	info.Bclr_used			= 0;
	info.BClr_important		= 0;
	fwrite(&info, sizeof(info), 1, image);

	for (int y = 0; y < GetH(); ++y) {
		const uint8_t* line = GetLine(info.Bheight - 1 - y);
		uint32_t padding = pitch - GetW() * pixelsize;
		for (int x = 0; x < GetW(); ++x, line += 4) {
			fwrite(line, pixelsize, 1, image);
		}
		for (int i = 0; i < static_cast<int>(padding); ++i)
			fputc(0, image);
	}
	fclose(image);
	return true;
}



uint32_t  Bitmap::SampleBilinear(float x, float y) const {
	int32_t fx = static_cast<int32_t>(x * 0x10000);
	int32_t fy = static_cast<int32_t>(y * 0x10000);
	int32_t x1 = Between(0, m_w - 1, fx >> 16);
	int32_t y1 = Between(0, m_h - 1, fy >> 16);
	int32_t x2 = Between(0, m_w - 1, x1 + 1);
	int32_t y2 = Between(0, m_h - 1, y1 + 1);
	int32_t dx = (fx >> 8) & 0xff;
	int32_t dy = (fy >> 8) & 0xff;

	if (m_w <= 0 || m_h <= 0) return 0;
	uint32_t c00 = GetPixel(x1, y1);
	uint32_t c01 = GetPixel(x2, y1);
	uint32_t c10 = GetPixel(x1, y2);
	uint32_t c11 = GetPixel(x2, y2);

	return BilinearInterp(c00, c01, c10, c11, dx, dy);
}

vec4f Bitmap::Sample2D(float u, float v) const {
	uint32_t res = SampleBilinear(u * m_w + 0.5f, v * m_h + .5f);
	return vector_from_color(res);
}

vec4f Bitmap::Sample2D(const vec2f& uv) const {
	return Sample2D(uv.x, uv.y);
}

void Bitmap::SetPixel(int x, int y, const vec4f& color) {
	SetPixel(x, y, vector_to_color(color));
}
void Bitmap::FlipVertical() {
	uint8_t* buffer = new uint8_t[m_pitch];
	for (int i = 0, j = m_h - 1; i < j; ++i, --j) {
		memcpy(buffer, GetLine(i), m_pitch);
		memcpy(GetLine(i), GetLine(j), m_pitch);
		memcpy(GetLine(j), buffer, m_pitch);
	}
	delete[] buffer;
}

void Bitmap::FlipHorizontal() {
	for (int y = 0; y < m_h; ++y) {
		for (int i = 0, j = m_w - 1; i < j; ++i, --j) {
			uint32_t c1 = GetPixel(i, y);
			uint32_t c2 = GetPixel(j, y);
			SetPixel(i, y, c2);
			SetPixel(j, y, c1);
		}
	}
}

uint32_t Bitmap::BilinearInterp(uint32_t tl, uint32_t tr,
	uint32_t bl, uint32_t br, int32_t distx, int32_t disty) {
	uint32_t f, r;
	int32_t distxy = distx * disty;
	int32_t distxiy = (distx << 8) - distxy;  
	int32_t distixy = (disty << 8) - distxy;  
	int32_t distixiy = 256 * 256 - (disty << 8) - (distx << 8) + distxy;
	r = (tl & 0x000000ff) * distixiy + (tr & 0x000000ff) * distxiy
		+ (bl & 0x000000ff) * distixy + (br & 0x000000ff) * distxy;
	f = (tl & 0x0000ff00) * distixiy + (tr & 0x0000ff00) * distxiy
		+ (bl & 0x0000ff00) * distixy + (br & 0x0000ff00) * distxy;
	r |= f & 0xff000000;
	tl >>= 16; tr >>= 16; bl >>= 16; br >>= 16; r >>= 16;
	f = (tl & 0x000000ff) * distixiy + (tr & 0x000000ff) * distxiy
		+ (bl & 0x000000ff) * distixy + (br & 0x000000ff) * distxy;
	r |= f & 0x00ff0000;
	f = (tl & 0x0000ff00) * distixiy + (tr & 0x0000ff00) * distxiy
		+ (bl & 0x0000ff00) * distixy + (br & 0x0000ff00) * distxy;
	r |= f & 0xff000000;
	return r;
}
