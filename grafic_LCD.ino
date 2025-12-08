#include <TFT_eSPI.h>
#define FOR(i,n) for(int i = 0;i < n;i++)
#define pai 3.141592
#define rgb tft.color565

TFT_eSPI tft = TFT_eSPI(); 
void setup() {
  tft.begin();
  tft.setRotation(0);

  tft.fillScreen(TFT_BLACK); 
}

void DrawPixel(int x,int y,uint16_t color){
  tft.drawPixel(x+tft.width()/2, -y+tft.height()/2, color);
}

void DrawLine(int x0,int y0,int x1,int y1,uint16_t color){
  tft.drawLine(x0+tft.width()/2, -y0+tft.height()/2, x1+tft.width()/2, -y1+tft.height()/2, color);
}
void DrawTriangle(int x0,int y0,int x1,int y1,int x2,int y2,uint16_t color){
  x0 = x0+tft.width()/2;
  x1 = x1+tft.width()/2;
  x2 = x2+tft.width()/2;
  y0 = -y0+tft.height()/2;
  y1 = -y1+tft.height()/2;
  y2 = -y2+tft.height()/2;
  tft.drawTriangle(x0,y0,x1,y1,x2,y2, color);
}
void FillTriangle(int x0,int y0,int x1,int y1,int x2,int y2,uint16_t color){
  x0 = x0+tft.width()/2;
  x1 = x1+tft.width()/2;
  x2 = x2+tft.width()/2;
  y0 = -y0+tft.height()/2;
  y1 = -y1+tft.height()/2;
  y2 = -y2+tft.height()/2;
  tft.fillTriangle(x0,y0,x1,y1,x2,y2, color);
}

void DrawFixedGUI(){
  int lineSpacing = 5;
  uint16_t scaleColor = TFT_GREEN;

  DrawLine(-5,0,5,0,scaleColor);//中心の水平線

  //目盛り生成
  FOR(i,10){//上半分
    if(i % 5 ==0){
      DrawLine(-15,lineSpacing*i,15,lineSpacing*i,scaleColor);
    }else{
      DrawLine(-10,lineSpacing*i,10,lineSpacing*i,scaleColor);
    }
  }
  FOR(i,10){//下半分
    if(i % 5 ==0){
      DrawLine(-15,lineSpacing*(-i),15,lineSpacing*(-i),scaleColor);
    }else{
      DrawLine(-10,lineSpacing*(-i),10,lineSpacing*(-i),scaleColor);
    }
  }

  //進行方向を示す三角形
  uint16_t dir_tri_color = TFT_WHITE;
  uint16_t dir_tri_outline_clr = TFT_GREEN;
  int ax = 13;
  int bx = ax+48;
  int cx = bx;
  int dx = cx - 19;
  int ay = 0;
  int by = ay-15;
  int cy = by-7;
  int dy = cy ;
  DrawTriangle(ax,ay,bx,by,dx,dy,dir_tri_color);
  DrawTriangle(bx,by,cx,cy,dx,dy,dir_tri_color);

  DrawTriangle(-ax,ay,-bx,by,-dx,dy,dir_tri_color);
  DrawTriangle(-bx,by,-cx,cy,-dx,dy,dir_tri_color);

  //画面の横真っ二つの線がわかるようにある三角形
  int ex = 100;
  int fx = ex+25;
  int ey = 0;
  int fy = 10;
  DrawTriangle(ex,ey,fx,fy,fx,-fy,dir_tri_color);
  DrawTriangle(-ex,ey,-fx,fy,-fx,-fy,dir_tri_color);

  //ロール目盛りの円弧
  for(int i = -100;i < 101;i++){
    DrawPixel(i,(int)sqrt(225.0 * 225.0 - (double)i * i);,dir_tri_outline_clr);
  }
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