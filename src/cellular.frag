precision highp float;                 // 浮動小数点精度を高く設定（演算誤差を低減）

uniform float uTime;                  // 経過時間
uniform vec2  u_resolution;            // 画面解像度（幅, 高さ）
uniform float num;                     // 球面上の点の総数
uniform float an;                      // カメラの回転角
out vec4 fragColor;                    // 出力するピクセルの最終 RGBA 値

// Inigo Quilez - Golden spiral sphere mapping demo

// p: 球面上の法線ベクトル（単位ベクトル）
// 戻り値: x = 最も近い黄金螺旋点のインデックス, y = その点との距離
vec2 inverseSF( vec3 p ) { 
    const float kTau = 6.28318530718;       // 2π
    const float kPhi = (1.0 + sqrt(5.0)) / 2.0; // 黄金比 Φ

    // 層（k）を推定
    // 高さ p.z から黄金螺旋ライン上のどの層にいるかを対数式で近似
    float k  = max(
        2.0, 
        floor(
            log2(num * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z))
            / log2(kPhi + 1.0)
        )
    );
    
    // フィボナッチ近似を使って F(k), F(k+1) を推定
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));   // (F(k), F(k+1))

    // 黄金螺旋グリッドの係数を生成
    vec2  ka = 2.0 * F / num;                       // 緯度方向の係数
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0)); // 経度方向の係数

    // 逆行列により点 p に最も近い格子番号 c を推定
    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x)
             / (ka.y * kb.x - ka.x * kb.y);

    vec2 c = floor(iB * vec2(
        atan(p.y, p.x),                  // 経度 φ
        p.z - 1.0 + 1.0 / num            // 緯度方向補正
    ));

    // 近傍の 4 点をチェックして最短距離の点を探す
    float d = 8.0;                       // 最短距離（初期値）
    float j = 0.0;                       // 最も近い点の index

    for (int s = 0; s < 4; s++) {
        vec2  uv = vec2(s & 1, s >> 1);  // (0,0),(1,0),(0,1),(1,1)
        float id = clamp(dot(F, uv + c), 0.0, num - 1.0);

        float phi = kTau * fract(id * kPhi);       // 経度
        float cosTheta = 1.0 - (2.0 * id + 1.0) / num; // 緯度 cosθ
        float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

        vec3 q = vec3(
            cos(phi) * sinTheta,
            sin(phi) * sinTheta,
            cosTheta
        );

        float tmp = dot(q - p, q - p);  // 距離の二乗
        if (tmp < d) {
            d = tmp;
            j = id;
        }
    }

    return vec2(j, sqrt(d));            // 点 index と距離
}

// 0~1 の乱数を返す簡易 hash 関数
float hash1(float n) { 
    return fract(sin(n) * 158.5453123);
}

vec2 hash22(vec2 p){
    float n = dot(p, vec2(127.1, 311.7));
    float m = dot(p, vec2(269.5, 183.3));
    return fract(sin(vec2(n, m)) * 43758.5453);
}

float fdist(vec2 p){
    vec2 n = floor(p + 0.5);
    float dist = sqrt(2.0);
    for(float j = -2.0; j <= 2.0; j++){
        for(float i = -2.0; i <= 2.0; i++){
            vec2 glid = n + vec2(i, j);
            vec2 jitter = sin(uTime) * (hash22(glid) - 0.5);
            dist = min(dist, distance(glid + jitter, p));
        }
    }
    return dist;
}

//---------------------------------------------------------------
// main
//---------------------------------------------------------------
void main() {

    // ピクセル座標を正規化（中心を (0,0) に）
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 p = (-u_resolution.xy + 2.0 * fragCoord.xy) / u_resolution.y;

    // カメラ設定
    vec3 ro = vec3(2.5 * cos(0.4*uTime), 1.0, 2.5 * sin(0.4*uTime));
    vec3 ta = vec3(0.0, 1.0, 0.0); // 注視点。カメラはここに向かっている

    vec3 ww = normalize(ta - ro); //カメラの前方向(正規化したta-ro)
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0))); //右方向ベクトル(wwとy軸との外積)
    vec3 vv = normalize(cross(uu, ww));  //上方向ベクトル(uuとwwの外積)
    //スクリーン上の座標p.xとp.yをカメラの右・上ベクトルに乗せ、前方向に1.5*wwを加えることで視野を作る
    //1.5は視野に相当する値
    vec3 rd = normalize(p.x * uu + p.y * vv + 1.5 * ww); 

    // 初期化
    vec3 sc = vec3(0.0, 1.0, 0.0); // 球の中心
    vec3 col = vec3(1.0);
    float tmin = 10000.0;         // 最小の交差距離(レイの最近接交差のt値(レイがどれだけ進んだか))
    //大きい初期値を設定しているのは、レイが何も当たらない場合を示すため
    vec3  nor = vec3(0.0);        //法線
    float occ = 1.0;              //影の係数
    vec3  pos = vec3(0.0);        //交差点


    // 地面 y=0 との交差判定
    float h = (0.0 - ro.y) / rd.y;
    if (h > 0.0) {
        tmin = h; //最初の交差距離
        nor = vec3(0.0, 1.0, 0.0); //地面の法線(上向き)
        pos = ro + h * rd; //交点座標の計算
        vec3 di = sc - pos; //地面上の交点から球中心までのベクトル(球との距離や影処理に使用)
        float l = length(di); //diの長さ
        occ = 1.0 - dot(nor, di / l) / (l * l);    // 影の処理
        col = vec3(1.0);
    }

    // 球との交差判定（2次方程式）
    vec3  ce = ro - sc;
    float b = dot(rd, ce);
    float c = dot(ce, ce) - 1.0;   
    h = b * b - c;                 //判別式

    if (h > 0.0) {
        h = -b - sqrt(h);
        if (h < tmin) {
            tmin = h;
            nor = normalize(ro + h * rd - sc); // 法線
            occ = 0.5 + 0.5 * nor.y;           // yが大きい(上向きほど明るい
        }

        //球面の法線をinverseSFに渡し、その方向上に最も近い螺旋状の点のインデックスとその点までの距離を返す
        vec2 fi = inverseSF(nor);


        float noise = fdist(nor.xy * 5.0 + fi.x * 13.7);
        float phase = uTime + fi.x * 0.2;
        float twinkle = 0.5 + 0.5 * sin(phase);

        float glow = exp(-40.0 * fi.y * fi.y);

        float rnd = hash1(fi.x * 13.0);
        vec3 starColor = mix(vec3(0.8,0.9,1.0), vec3(0.7,0.8,1.0), rnd);

        col = starColor * (0.4 + glow * twinkle) * (0.5 + 0.5 * noise);
    }
    
    // 背景フェード
    if (tmin < 100.0) {
        pos = ro + tmin * rd;
        col *= occ;
        col = mix(col, vec3(1.0), 1.0 - exp(-0.003 * tmin * tmin));
    }

    // ガンマ補正
    col = sqrt(col);
    fragColor = vec4(col, 1.0);
}
