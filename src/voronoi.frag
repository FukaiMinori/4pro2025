precision highp float;

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

float fdist21(vec2 p, float freq){
    vec2 n = floor(p + 0.5); 
    float dist = sqrt(2.0);
    
    for(float j = 0.0; j <= 2.0; j++){
        vec2 glid;
        glid.y = n.y + sign(mod(j, 2.0) - 0.5) * ceil(j * 0.5);
        if(abs(glid.y - p.y) - 0.5 > dist){
            continue;
        }
        for(float i = -1.0; i <= 1.0; i++){
            glid.x = n.x + i;
            
            // ★ポイント: X座標を周波数でループさせる
            // これにより uv.x=0.0 と uv.x=1.0 が同じハッシュ値を参照し、繋ぎ目が消える
            vec2 wrappedGlid = vec2(mod(glid.x, freq), glid.y);
            
            vec2 jitter = hash22(wrappedGlid) - 0.5;
            dist = min(dist, length(glid + jitter - p));
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

    //カメラ位置
    vec3 ro = vec3(2.5 * cos(an), 1.0, 2.5 * sin(an));
    vec3 ta = vec3(0.0, 1.0, 0.0); // 注視点

    vec3 ww = normalize(ta - ro);
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0)));
    vec3 vv = normalize(cross(uu, ww));
    vec3 rd = normalize(p.x * uu + p.y * vv + 1.5 * ww); 

    // 初期化
    vec3 sc = vec3(0.0, 1.0, 0.0); // 球の中心
    vec3 col = vec3(1.0);
    float tmin = 10000.0;
    vec3  nor = vec3(0.0);        //法線
    float occ = 1.0;
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
            nor = normalize(ro + h * rd - sc); 
            occ = 0.5 + 0.5 * nor.y;
        }

        vec2 fi = inverseSF(nor);

        const float TAU = 6.28318530718;
        const float PI = 3.14159265359;
        
        // 経度を -PI〜PI から 0〜1 に変換
        float phi = atan(nor.y, nor.x); 
        float theta = acos(nor.z);
        vec2 uv = vec2(phi / TAU + 0.5, theta / PI);

        // ★ポイント: freq を整数にキャストする（1周で模様を完結させるため）
        float freq = floor(sqrt(max(num, 1.0)));

        float dNorm = fi.y * sqrt(num);
        float fibMask = 1.0 - smoothstep(0.2, 1.0, dNorm);
        float fibRand = hash1(fi.x * 7.123);

        vec2 fibOffset = 0.05 * fibMask * vec2(
            sin(fibRand * 6.2831 + uTime * 0.2),
            cos(fibRand * 6.2831 + uTime * 0.2)
        );

        // ★修正した fdist21 を呼び出す（引数に freq を追加）
        float n = fdist21(uv * freq + fibOffset, freq);

        float gray = 1.0 - n;
        gray = pow(gray, 0.7);
        col = vec3(gray);
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
