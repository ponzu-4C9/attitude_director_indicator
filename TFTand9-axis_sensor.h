#include <TFT_eSPI.h>
#include <Wire.h>
#include <ICM20948_WE.h>
#define ICM20948_ADDR 0x68
#define rgb tft.color565

class TFTand9axis_sensor{

  private:
  TFT_eSPI tft;
  ICM20948_WE myIMU;

  typedef struct {
    double x;
    double y;
    double z;
    double w;
  } Vector;

  typedef struct {
    Vector a;
    Vector b;
  } Line;

  typedef struct {
    double x;
    double y;
  } Point;


  void applyMatrix(double a, double b, double c, double d, double *x, double *y) {

    // 計算結果を一時的に保持する変数
    double newX;
    double newY;

    // 新しい X 成分 = a*x + b*y
    newX = a * *x + b * *y;

    // 新しい Y 成分 = c*x + d*y
    newY = c * *x + d * *y;

    // 結果を出力ポインタに書き込み
    *x = newX;
    *y = newY;
  }

  //描画関数を数学的な座標にオフセット
  void DrawPixel(int x, int y, uint16_t color) {
    tft.drawPixel(x + tft.width() / 2, -y + tft.height() / 2, color);
  }

  void DrawLine(int x0, int y0, int x1, int y1, uint16_t color) {
    tft.drawLine(x0 + tft.width() / 2, -y0 + tft.height() / 2, x1 + tft.width() / 2, -y1 + tft.height() / 2, color);
  }
  void DrawTriangle(int x0, int y0, int x1, int y1, int x2, int y2, uint16_t color) {
    x0 = x0 + tft.width() / 2;
    x1 = x1 + tft.width() / 2;
    x2 = x2 + tft.width() / 2;
    y0 = -y0 + tft.height() / 2;
    y1 = -y1 + tft.height() / 2;
    y2 = -y2 + tft.height() / 2;
    tft.drawTriangle(x0, y0, x1, y1, x2, y2, color);
  }
  void FillTriangle(int x0, int y0, int x1, int y1, int x2, int y2, uint16_t color) {
    x0 = x0 + tft.width() / 2;
    x1 = x1 + tft.width() / 2;
    x2 = x2 + tft.width() / 2;
    y0 = -y0 + tft.height() / 2;
    y1 = -y1 + tft.height() / 2;
    y2 = -y2 + tft.height() / 2;
    tft.fillTriangle(x0, y0, x1, y1, x2, y2, color);
  }

  //ロール角を表す三角形の初期化
  Point topTri[3];
  const double gy_dis = 200;  //原点からの円弧の距離
  void topTriinitializer(){
    const int size = 20;
    topTri[0].x = 0;topTri[0].y = gy_dis;
    topTri[1].x = size/2;topTri[1].y = gy_dis-size;
    topTri[2].x = -size/2;topTri[2].y = gy_dis-size;
  }

  //目盛り等の生成
  void DrawFixedGUI(int lineSpacing = 5) {

    uint16_t scaleColor = TFT_WHITE;

    DrawLine(-5, 0, 5, 0, scaleColor);  //中心の水平線

    //目盛り生成
    for (int i = 1; i <= 10; i++) {  //上半分
      if (i % 5 == 0) {
        DrawLine(-15, lineSpacing * i, 15, lineSpacing * i, scaleColor);
      } else {
        DrawLine(-10, lineSpacing * i, 10, lineSpacing * i, scaleColor);
      }
    }
    for (int i = 1; i <= 10; i++) {  //下半分
      if (i % 5 == 0) {
        DrawLine(-15, lineSpacing * (-i), 15, lineSpacing * (-i), scaleColor);
      } else {
        DrawLine(-10, lineSpacing * (-i), 10, lineSpacing * (-i), scaleColor);
      }
    }

    //進行方向を示す三角形
    uint16_t dir_tri_color = TFT_WHITE;
    uint16_t dir_tri_outline_clr = TFT_WHITE;
    int ax = 13;
    int bx = ax + 48;
    int cx = bx;
    int dx = cx - 19;
    int ay = 0;
    int by = ay - 15;
    int cy = by - 7;
    int dy = cy;
    DrawTriangle(ax, ay, bx, by, dx, dy, dir_tri_outline_clr);
    DrawTriangle(bx, by, cx, cy, dx, dy, dir_tri_outline_clr);

    DrawTriangle(-ax, ay, -bx, by, -dx, dy, dir_tri_outline_clr);
    DrawTriangle(-bx, by, -cx, cy, -dx, dy, dir_tri_outline_clr);

    //画面の横真っ二つの線がわかるようにある三角形
    uint16_t baseTriangleColor = TFT_WHITE;
    int ex = 100;
    int fx = ex + 25;
    int ey = 0;
    int fy = 10;
    DrawTriangle(ex, ey, fx, fy, fx, -fy, baseTriangleColor);
    DrawTriangle(-ex, ey, -fx, fy, -fx, -fy, baseTriangleColor);


    //ロール目盛りの生成
    uint16_t rollScaleColor = TFT_WHITE;
    int scaleLength = 10;
    const double scaleIntervalAngle = 5 * (PI / 180);  //度に(pai/180)を掛けるとラジアンに変換できる。
    int tickCount = 7;


    double gx = 0;
    double gy = gy_dis;
    double hx = 0;
    double hy = gy + scaleLength;
    double startAngle = -scaleIntervalAngle * (tickCount + 1);
    applyMatrix(cos(startAngle), -sin(startAngle), sin(startAngle), cos(startAngle), &gx, &gy);
    applyMatrix(cos(startAngle), -sin(startAngle), sin(startAngle), cos(startAngle), &hx, &hy);
    FOR(i, tickCount * 2 + 1) {
      applyMatrix(cos(scaleIntervalAngle), -sin(scaleIntervalAngle), sin(scaleIntervalAngle), cos(scaleIntervalAngle), &gx, &gy);
      applyMatrix(cos(scaleIntervalAngle), -sin(scaleIntervalAngle), sin(scaleIntervalAngle), cos(scaleIntervalAngle), &hx, &hy);
      DrawLine(gx, gy, hx, hy, rollScaleColor);
    }

    //ロール目盛りの円弧
    int dar = gy_dis * sin(scaleIntervalAngle * tickCount);  //一番端っこの目盛りの付け根のx座標
    for (int i = -dar; i < dar; i++) {
      DrawPixel(i, (int)sqrt((gy_dis * gy_dis) - (double)i * i), rollScaleColor);
    }
  }



  #define SIZE 4
  void multiply_matrix(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double C[SIZE][SIZE]) {
    // 行列の積 C[i][j] = Sum(A[i][k] * B[k][j])
    for (int i = 0; i < SIZE; i++) {      // 結果Cの行 (i)
      for (int j = 0; j < SIZE; j++) {    // 結果Cの列 (j)
        C[i][j] = 0.0f;                   // 初期化
        for (int k = 0; k < SIZE; k++) {  // 共通のインデックス (k)
          C[i][j] += A[i][k] * B[k][j];
        }
      }
    }
  }

  void applyMatrix4(const double M[SIZE][SIZE], Vector *V) {
    Vector temp;
    temp.x = M[0][0] * V->x + M[0][1] * V->y + M[0][2] * V->z + M[0][3] * V->w;
    temp.y = M[1][0] * V->x + M[1][1] * V->y + M[1][2] * V->z + M[1][3] * V->w;
    temp.z = M[2][0] * V->x + M[2][1] * V->y + M[2][2] * V->z + M[2][3] * V->w;
    temp.w = M[3][0] * V->x + M[3][1] * V->y + M[3][2] * V->z + M[3][3] * V->w;
    V->x = temp.x;
    V->y = temp.y;
    V->z = temp.z;
    V->w = temp.w;
  }

  double z_max = 1000;
  const int h_num = 10;
  const double hrs = z_max / h_num;  //横線間隔
  const double vrs = hrs/10;//縦線間隔
  const int v_num = 11;//縦線は奇数個にしてそうじゃないと中心にこない
  const double ground_width = 200;

  Line hline[h_num];
  void hlineinitializer() {  //これは直線だからとりあえず長さを1の線分として定義
    for (int i = 0; i < h_num; i++) {
      hline[i].a.x = -ground_width/2;
      hline[i].b.x = ground_width/2;
      hline[i].a.y = hline[i].b.y = 0;
      hline[i].a.z = hline[i].b.z = (i+1) * hrs;
      hline[i].a.w = hline[i].b.w = 1;
    }
  }
  Line vline[v_num];  //始点z_max,終点0の線分
  void vlineinitializer() {
    for (int i = 0; i < v_num; i++) {
      vline[i].a.x = vline[i].b.x = (i-(v_num/2)) * vrs;
      vline[i].a.y = vline[i].b.y = 0;
      vline[i].a.z = z_max;
      vline[i].b.z = 10;
      vline[i].a.w = vline[i].b.w = 1;
    }
  }
  double T[SIZE][SIZE] = {
    { 1, 0, 0, 0 },
    { 0, 1, 0, 0 },
    { 0, 0, 1, 0 },
    { 0, 0, 0, 1 }
  } ;
  double Rx[SIZE][SIZE] = {
    { 1, 0, 0, 0 },
    { 0, 1, 0, 0 },
    { 0, 0, 1, 0 },
    { 0, 0, 0, 1 }
  } ;
  double Rz[SIZE][SIZE] = {
    { 1, 0, 0, 0 },
    { 0, 1, 0, 0 },
    { 0, 0, 1, 0 },
    { 0, 0, 0, 1 }
  } ;

  double Rz_tri[2][2];

  double M_temp[SIZE][SIZE];
  double M[SIZE][SIZE];

  Line hlineTemp[h_num];
  Line vlineTemp[v_num];

  int hxy[h_num][4];
  int vxy[v_num][4];

  public:
  void begin(){
    tft = TFT_eSPI();
    myIMU = ICM20948_WE(ICM20948_ADDR);
    Wire.begin();
    if (!myIMU.init()) {
      Serial.println("ICM20948 does not respond");
    } else {
      Serial.println("ICM20948 is connected");
    }
    myIMU.autoOffsets();
    myIMU.setAccRange(ICM20948_ACC_RANGE_2G);
    myIMU.setAccDLPF(ICM20948_DLPF_6);
    myIMU.setAccSampleRateDivider(10);

    tft.begin();
    tft.setRotation(0);
    tft.fillScreen(TFT_BLACK);
  }
  float getPitchAndRoll(float *p,float *r){
    myIMU.readSensor();
    *p = myIMU.getPitch();
    *r = myIMU.getRoll();
  }
  void updata(double h/*高度(px)*/,){

    float pitch;
    float roll;
    getPitchAndRoll(&pitch,&roll);

    double Spacing = 10;

    double f = 1000;//3Dグラフィックの焦点

    double fai = pitch * (PI / 180);
    double theta = roll * (PI / 180);

    Serial.printf("pitch%.15f, roll%.15f\n",pitch,roll);


    Rz_tri[0][0] = cos(-theta);Rz_tri[0][1] = -sin(-theta);
    Rz_tri[1][0] = sin(-theta);Rz_tri[1][1] = cos(-theta);
    DrawTriangle(topTri[0].x,topTri[0].y,topTri[1].x,topTri[1].y,topTri[2].x,topTri[2].y,TFT_BLACK);
    topTriinitializer();
    for(int i = 0;i < 3;i++){
      applyMatrix(Rz_tri[0][0],Rz_tri[0][1],Rz_tri[1][0],Rz_tri[1][1],&topTri[i].x,&topTri[i].y);
    }
    DrawTriangle(topTri[0].x,topTri[0].y,topTri[1].x,topTri[1].y,topTri[2].x,topTri[2].y,TFT_GREEN);


    double Spfai = Spacing*pitch;
    double numerator = -(Spfai * z_max + h * f * cos(theta));
    double denominator = z_max * f * cos(theta) - Spfai * h;
    double phi = atan2(numerator, denominator);

    T[1][3] = -h;  //下にh文下げる

    Rx[1][1] = cos(-phi);Rx[1][2] = sin(-phi);
    Rx[2][1] = -sin(-phi);Rx[2][2] = cos(-phi);

    Rz[0][0] = cos(-theta);Rz[1][0] = sin(-theta);
    Rz[0][1] = -sin(-theta);Rz[1][1] = cos(-theta);


    multiply_matrix(Rz,Rx,M_temp);
    multiply_matrix(M_temp,T,M);

    for (int i = 0; i < h_num; i++) {
      DrawLine(hxy[i][0], hxy[i][1], hxy[i][2], hxy[i][3], TFT_BLACK);
    }
    for (int i = 0; i < v_num; i++) {
      DrawLine(vxy[i][0], vxy[i][1], vxy[i][2], vxy[i][3], TFT_BLACK);
    }
    
    hlineinitializer();
    vlineinitializer();

    for (int i = 0; i < h_num; i++) {
      applyMatrix4(M,&(hline[i].a));
      applyMatrix4(M,&(hline[i].b));
      if (hline[i].a.z > 1.0 && hline[i].b.z > 1.0) { // カメラより後ろ（z<=0）を描画しないガード
        hxy[i][0] = (f * hline[i].a.x / hline[i].a.z);
        hxy[i][1] = (f * hline[i].a.y / hline[i].a.z);
        hxy[i][2] = (f * hline[i].b.x / hline[i].b.z);
        hxy[i][3] = (f * hline[i].b.y / hline[i].b.z);
      }else {
        // 【重要】Zが手前すぎたら座標をリセットして、描画（消去）させない
        hxy[i][0] = hxy[i][1] = hxy[i][2] = hxy[i][3] = 0;
      }
    }

    for (int i = 0; i < v_num; i++) {
      applyMatrix4(M,&(vline[i].a));
      applyMatrix4(M,&(vline[i].b));
      if (vline[i].a.z > 0 && vline[i].b.z > 0) { // カメラより後ろ（z<=0）を描画しないガード
        vxy[i][0] = (f * vline[i].a.x / vline[i].a.z);
        vxy[i][1] = (f * vline[i].a.y / vline[i].a.z);
        vxy[i][2] = (f * vline[i].b.x / vline[i].b.z);
        vxy[i][3] = (f * vline[i].b.y / vline[i].b.z);
      }
    }

    int haaa[4];
    for (int i = 0; i < h_num; i++) {
      DrawLine(hxy[i][0], hxy[i][1], hxy[i][2], hxy[i][3], TFT_GREEN);
      if(i == h_num-1){
        haaa[0]=hxy[i][0];
        haaa[1]=hxy[i][1];
        haaa[2]=hxy[i][2];
        haaa[3]=hxy[i][3];
      }
    }
    double a = (double)(haaa[3] - haaa[1]) / (double)(haaa[2] - haaa[0]);
    double b = haaa[1] - a * haaa[0];

    for(int j = 0;j < 4;j++){
    Serial.printf("%d|", haaa[j]);
    }
    Serial.printf("\ny = %f * x + %f\n", a, b);


    for (int i = 0; i < v_num; i++) {
      DrawLine(vxy[i][0], vxy[i][1], vxy[i][2], vxy[i][3], TFT_GREEN);
    }

    DrawFixedGUI(Spacing);
  }
}