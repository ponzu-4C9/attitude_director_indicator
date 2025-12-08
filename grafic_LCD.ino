#include <TFT_eSPI.h>
#define FOR(i,n) for(int i = 0;i < n;i++)
#define pai 3.141592

uint16_t rgb(uint8_t r, uint8_t g, uint8_t b)
{
  // 1. R成分を5ビットに変換し、11ビット左にシフト (R_4 R_3 ... R_0 | 0000 0000 000)
  // 255を31に圧縮するため、右に3ビットシフト (R_8bit >> 3)
  uint16_t r5 = (r >> 3) << 11;

  // 2. G成分を6ビットに変換し、5ビット左にシフト (0000 0 | G_5 G_4 ... G_0 | 00000)
  // 255を63に圧縮するため、右に2ビットシフト (G_8bit >> 2)
  uint16_t g6 = (g >> 2) << 5;

  // 3. B成分を5ビットに変換 (0000 0000 000 | B_4 B_3 ... B_0)
  // 255を31に圧縮するため、右に3ビットシフト (B_8bit >> 3)
  uint16_t b5 = b >> 3;

  // 4. 全ての成分をビットOR演算で結合
  return r5 | g6 | b5;
}

void DrawLine(int x0,int y0,int x1,int y1,uint16_t color){
  int cw = tft.width()/2;
  int ch = tft.height()/2;
  tft.drawLine(x0+cw, -y0+ch, x1+cw, -y1+ch, color);
}

TFT_eSPI tft = TFT_eSPI(); 


void setup() {
  tft.begin();
  tft.setRotation(0);
  cw = tft.width()/2;
  ch = tft.height()/2;

  tft.fillScreen(TFT_BLACK); 
}



double theta = 0;
int y0a = 0;
int y1a = 0;
void loop() {
  DrawLine(-160,y0a,160,y1a, TFT_BLACK);
  double fai = 0;
  double k = 0;
  theta += 0.001;

  y0a = tan(theta)*(-160)+(-k*fai);
  y1a = tan(theta)*(160)+(-k*fai);

  DrawLine(-160,y0a,160,y1a, rgb(105,212,16));
  
}