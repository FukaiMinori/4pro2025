
precision highp float;//浮動小数点の精度を「高」に設定(グラフィックスの計算誤差を減らすため)

uniform float u_time;//経過時間
uniform vec2  u_resolution;//画面サイズ(幅・高さ) Shadertoy互換名に統一
uniform float radius;
uniform float num; //球面上に配置する総点数(フィボナッチ格子の点数)
uniform float an; //カメラ位置の回転
out vec4 fragColor;//各ピクセルの最終的な出力色

// ---------------------------------------------------------------
// Inigo Quilez - https://iquilezles.org/
// Golden spiral sphere mapping demo
// ---------------------------------------------------------------

//関数：球面上の点p(x,y,z)が黄金比を用いた点列の中でどのインデックス(番号)の点に最も近いか
vec2 inverseSF( vec3 p ) { 
    const float kTau = 6.28318530718; //2π
    const float kPhi = (1.0 + sqrt(5.0)) / 2.0; //黄金比Φ
    //const float num = 130.0; //球面上の点の総数

    //球面上の点の「層」を求める。高さp.zに応じて螺旋のどの位置にいるかを推定
    float k  = max(2.0, floor(log2(num * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z)) / log2(kPhi + 1.0)));
    //フィボナッチ数近似式を使ってn番目とその次のフィボナッチ数を推定
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi)); //Fk=n番目, Fk*kPhi=n+1番目

    //黄金螺旋上の格子を球面にマッピングするための係数を生成
    vec2  ka = 2.0 * F / num;
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));

    //逆行列を使ってpの座標に最も近い「格子点」を見つける
    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);
    vec2 c = floor(iB * vec2(atan(p.y, p.x), p.z - 1.0 + 1.0 / num));

    //最短距離を初期化
    float d = 8.0;
    float j = 0.0;
    for (int s = 0; s < 4; s++) { //周囲4点(近傍を格子)を調べて、最も近い点を選ぶ
        //ビット演算で(0,0),(1,0),(0,1),(1,1) の近傍をループ
        vec2  uv = vec2(s & 1, s >> 1);
        float id = clamp(dot(F, uv + c), 0.0, num - 1.0);

        //球面座標に変換
        float phi = kTau * fract(id * kPhi);
        float cosTheta = 1.0 - (2.0 * id + 1.0) / num;
        float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

        //実際の点pとqの距離を測り、最も近い点のIDを保存
        vec3 q = vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
        float tmp = dot(q - p, q - p);
        if (tmp < d) {
            d = tmp;
            j = id;
        }
    }
    return vec2(j, sqrt(d));//その点のインデックス(j)と距離(d)
}

//hash1()：乱数生成→簡易ノイズ関数。nから一様乱数(0~1)を作る
float hash1(float n) { return fract(sin(n) * 158.5453123); }

//描画
void main() {
    //ピクセル座標を正規化して、中心を(0,0)にする
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 p = (-u_resolution.xy + 2.0 * fragCoord.xy) / u_resolution.y;

    //カメラ設定
    //float an = 0.5 * u_time;
    //float an = -3.14159 / 2.0;
    vec3 ro = vec3(2.5 * cos(an), 1.0, 2.5 * sin(an)); //カメラ位置
    vec3 ta = vec3(0.0, 1.0, 0.0); //注視点
    vec3 ww = normalize(ta - ro);
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0)));
    vec3 vv = normalize(cross(uu, ww));
    vec3 rd = normalize(p.x * uu + p.y * vv + 1.5 * ww); //レイ(光線)の方向ベクトル

    //レイマーチング(光線追跡)
    vec3 sc = vec3(0.0, 1.0, 0.0); //球の中心を設定
    vec3 col = vec3(1.0); //色初期化

    //最近交差点情報を初期化
    float tmin = 10000.0;
    vec3  nor = vec3(0.0);
    float occ = 1.0;
    vec3  pos = vec3(0.0);

    //地面との交差判定
    //y=0の地面にレイが当たったときの処理。明るさや影(occlusion)も計算
    float h = (0.0 - ro.y) / rd.y;
    if (h > 0.0) {
        tmin = h;
        nor = vec3(0.0, 1.0, 0.0);
        pos = ro + h * rd;
        vec3 di = sc - pos;
        float l = length(di);
        occ = 1.0 - dot(nor, di / l) * 1.0 * 1.0 / (l * l);
        col = vec3(1.0);
    }

    //球との交差判定
    //球の方程式を解く(レイ交差判定)
    vec3  ce = ro - sc;
    float b = dot(rd, ce);
    float c = dot(ce, ce) - 1.0;
    h = b * b - c;
    //球とレイが交差したら、法線ベクトルを求める
    if (h > 0.0) {
        h = -b - sqrt(h);
        if (h < tmin) {
            tmin = h;
            nor = normalize(ro + h * rd - sc);
            occ = 0.5 + 0.5 * nor.y;
        }

        //fi.x=点のインデックス(どのスパイラル位置か) fi.y=近さ(距離)(点の位置に応じて色をつける)
        vec2 fi = inverseSF(nor);
        col = 0.5 + 0.5 * sin(hash1(fi.x * 13.0) * 3.0 + 1.0 + vec3(0.0, 1.0, 1.0));
        col *= smoothstep(0.001*num, 0.005/num, fi.y);
        col *= mix(1.0, 1.0 - smoothstep(0.12, 0.125, fi.y), smoothstep(-0.9, 0.1, sin(u_time)));
        col *= 1.0 + 0.1 * sin(100.0*radius * fi.y);
        col *= 1.5;
    }

    //背景・減衰処理
    //奥行き(tmin)に応じて色をフェードアウト
    if (tmin < 100.0) {
        pos = ro + tmin * rd;
        col *= occ;
        col = mix(col, vec3(1.0), 1.0 - exp(-0.003 * tmin * tmin));
    }

    //ガンマ補正と出力
    //色を人間の視覚に合わせて補正し、最終出力
    col = sqrt(col);
    fragColor = vec4(col, 1.0);
}
