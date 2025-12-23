precision highp float; // 浮動小数点精度を高く設定（演算誤差を低減）

uniform vec2  u_resolution; // 画面解像度（幅, 高さ）
uniform float num; // 球面上の点の総数
uniform float an; // カメラの回転角
uniform float embedSmall;
uniform float embedSmall2;
out vec4 fragColor; // 出力するピクセルの最終 RGBA 値

// Inigo Quilez - Golden spiral sphere mapping demo

// p: 球面上の法線ベクトル（単位ベクトル）
// 戻り値:x = 最も近い黄金螺旋点のインデックス, y = その点との距離
vec2 inverseSF( vec3 p ) { 
    const float kTau = 6.28318530718; // 2π
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
    vec2  ka = 2.0 * F / num;// 緯度方向の係数
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0)); // 経度方向の係数

    // 逆行列により点 p に最も近い格子番号 c を推定
    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x)
             / (ka.y * kb.x - ka.x * kb.y);

        vec2 c = floor(iB * vec2(
        atan(p.y, p.x), // 経度 φ
        p.z - 1.0 + 1.0 / num // 緯度方向補正
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

vec3 cauliflowerColor(float seed, float shade) {
    // ベース：白〜黄緑
    vec3 base = mix(
        vec3(0.95, 0.96, 0.90),   // 白
        vec3(0.8039, 0.9059, 0.7373),   // 黄緑
        shade
    );

    // 微妙なムラ
    float n = hash1(seed * 17.3);
    base += 0.05 * (n - 0.5);

    return clamp(base, 0.0, 1.0);
}

// インデックス id から単位球面上の3D位置 q を返す
vec3 fibonacciPoint(float id, float total) {
    const float kTau = 6.28318530718;
    const float kPhi = (1.0 + sqrt(5.0)) / 2.0;
    float phi = kTau * fract(id * kPhi);
    float cosTheta = 1.0 - (2.0 * id + 1.0) / total;
    float sinTheta = sqrt(max(0.0, 1.0 - cosTheta * cosTheta));
    return vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

bool intersectSphere(vec3 ro, vec3 rd, vec3 center, float radius, out float t) {
    vec3 oc = ro - center;
    float b = dot(rd, oc);
    float c = dot(oc, oc) - radius * radius;
    float D = b*b - c;
    if (D < 0.0) return false;
    float s = sqrt(max(D, 0.0)); // クランプ
    float t0 = -b - s;
    float t1 = -b + s;
    t = (t0 > 0.0) ? t0 : ((t1 > 0.0) ? t1 : -1.0);
    return t > 0.0;
}
// 交差比較用の相対EPS
float relEps(float t) {
    return max(1e-6, 1e-5 * max(t, 1.0));
}

// 近傍4候補の id を返す（inverseSF の c 推定を再利用）
void getNearest4Ids(vec3 p, out float ids[4]) {
    const float kTau = 6.28318530718;
    const float kPhi = (1.0 + sqrt(5.0)) / 2.0;

    float k = max(2.0,
        floor( log2(num * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z)) / log2(kPhi + 1.0) )
    );

    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi)); // (F(k), F(k+1))
    vec2 ka = 2.0 * F / num;
    vec2 kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));

    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);

    vec2 c = floor(iB * vec2( atan(p.y, p.x), p.z - 1.0 + 1.0 / num ));

    // 4近傍 (0,0),(1,0),(0,1),(1,1)
    for (int s = 0; s < 4; s++) {
        vec2 uv = vec2(s & 1, s >> 1);
        float id = clamp(dot(F, uv + c), 0.0, num - 1.0);
        ids[s] = id;
    }
}

//子球処理
bool hitSubSphere(
    vec3 ro, vec3 rd,
    vec3 center, float radius,
    float tBigOcc, bool bigBlocks,
    inout float tmin, inout vec3 pos, inout vec3 col, inout float occ
){
    float t;
    if (!intersectSphere(ro, rd, center, radius, t)) return false;
    if (t <= 0.0) return false;

    // 大球の遮蔽チェック
    if (bigBlocks && t >= tBigOcc - relEps(t)) return false;
    if (t >= tmin - relEps(tmin)) return false;

    // ヒット情報
    vec3 hitPos = ro + t * rd;
    vec3 hitNor = normalize(hitPos - center);
    vec2 fi     = inverseSF(hitNor);

    float shade = clamp(0.6 + 0.4 * fi.y, 0.0, 1.0);
    vec3  c     = cauliflowerColor(fi.x + 10.0, shade);

    // 子球は少し明るく
    c *= 1.05;

    float o = 0.5 + 0.5 * hitNor.y;

    // 更新
    tmin = t;
    pos  = hitPos;
    col  = c;
    occ  = mix(1.0, o, 0.35);

    return true;
}

// 任意：大球外突起のフィルタ（使う場合のみ呼ぶ）
bool protrudesOut(vec3 center, float r, vec3 bigCenter, float bigR) {
    return (length(center - bigCenter) + r) > (bigR + 1e-5);
}

void testBump(vec3 ro, vec3 rd, vec3 sc, vec3 q,
   inout float tmin, inout vec3 pos, inout vec3 col, inout float occ)
{
    float bumpR0    = 0.3;
    float subRScale = 0.35;

    float rS      = bumpR0;
    float dCenter = 1.0 + rS * (1.0 - embedSmall);
    vec3  sphC    = sc + q * dCenter;

    float subR    = bumpR0 * subRScale;

    float tBigOcc;
    bool bigBlocks = intersectSphere(ro, rd, sc, 1.0, tBigOcc);

    // --- 子球（初期方向） ---
    vec3 q2_init   = normalize(q);
    vec3 subC_init = sphC + q2_init * (rS + subR * (1.0 - embedSmall2));

    hitSubSphere(ro, rd, subC_init, subR, tBigOcc, bigBlocks,
                 tmin, pos, col, occ);

    // --- 親球 ---
    float tS;
    if (intersectSphere(ro, rd, sphC, rS, tS) && tS > 0.0 &&
        (!bigBlocks || tS < tBigOcc - relEps(tS)) &&
        tS < tmin - relEps(tmin))
    {
        vec3 hitPos = ro + tS * rd;
        vec3 hitNor = normalize(hitPos - sphC);
        vec2 fiHit  = inverseSF(hitNor);

        float shade = clamp(0.4 + 0.6 * fiHit.y, 0.0, 1.0);
        vec3 iterCol = cauliflowerColor(fiHit.x, shade);

        iterCol *= 0.99 + 0.20 * hitNor.y;
        float iterOcc = 0.2 + 0.2 * hitNor.y;

        tmin = tS;
        pos  = hitPos;
        col  = iterCol;
        occ  = mix(1.0, iterOcc, 0.35);

        // --- 精密化子球 ---
        vec3 q2_ref   = normalize(fibonacciPoint(fiHit.x, 1.0));
        vec3 subC_ref = sphC + q2_ref * (rS + subR * (1.0 - embedSmall2));

        hitSubSphere(ro, rd, subC_ref, subR, tBigOcc, bigBlocks,
                     tmin, pos, col, occ);
    }
}

// main: ピクセルごとの描画処理
void main() {
    // 正規化スクリーン座標
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 p = (-u_resolution.xy + 2.0 * fragCoord.xy) / u_resolution.y;

    // カメラ
    vec3 ro = vec3(2.5 * cos(an), 1.0, 2.5 * sin(an));
    vec3 ta = vec3(0.0, 1.0, 0.0);
    vec3 ww = normalize(ta - ro);
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0)));
    vec3 vv = normalize(cross(uu, ww));
    vec3 rd = normalize(p.x * uu + p.y * vv + 1.5 * ww);

    // 初期化
    vec3 sc = vec3(0.0, 1.0, 0.0);   // 大球中心
    vec3 col = vec3(0.0);            // 出力色
    float tmin = 1e20;                // 最前面距離
    float occ  = 1.0;                // 簡易オクルージョン
    const float EPS = 1e-6;
    vec3  pos = vec3(0.0);
    vec3  nor = vec3(0.0);

    // 大球（半径1）
    vec3 ce = ro - sc;
    float b = dot(rd, ce);
    float c = dot(ce, ce) - 1.0;

    vec2 fi = vec2(0.0);

    // 大球
    // 大球ヒットで fi・nor を決定
    float hBig;
    bool hitBig = false;
    if (intersectSphere(ro, rd, sc, 1.0, hBig)) {
        if (hBig > 0.0 && hBig + relEps(tmin) < tmin) {
            tmin = hBig;
            pos  = ro + hBig * rd;
            nor  = normalize(pos - sc);
            occ  = 0.5 + 0.5 * nor.y;
            col = vec3(0.92, 0.94, 0.88) * occ;
            fi   = inverseSF(nor);
            hitBig = true;
        }
    }

 // 背景側: hitBig == false のとき
if (!hitBig) {
    float tBigOcc;
    bool bigBlocks = intersectSphere(ro, rd, sc, 1.0, tBigOcc); // 大球の潜在的オクルージョン距離

    vec3 pClosest = ro + rd * max(-dot(rd, ro - sc), 0.0);
    vec3 probeDir = normalize(pClosest - sc);

    float ids[4];
    getNearest4Ids(probeDir, ids);
    for (int i=0; i<4; i++) {
        vec3 q = fibonacciPoint(ids[i], num);

        // 親・子の中心と半径を「事前に」作って外突起か確認
        float bumpR0    = 0.3;
        float subRScale = 0.35;
        float rS        = bumpR0;
        float dCenter   = 1.0 + rS * (1.0 - embedSmall);
        vec3  sphC      = sc + q * dCenter;
        float subR      = bumpR0 * subRScale;

        vec3 q2_init    = normalize(q);
        vec3 subC_init  = sphC + q2_init * (rS + subR * (1.0 - embedSmall2));

        // 外向き & 大球外へ突き出しているものだけ許可
        if (dot(q, probeDir) > 0.0 && protrudesOut(subC_init, subR, sc, 1.0)) {
            testBump(ro, rd, sc, q, tmin, pos, col, occ);
        }
    }
}

// 大球ヒット側: hitBig == true のとき
if (hitBig) {
    float ids[4];
    getNearest4Ids(nor, ids);
    for (int i=0; i<4; i++) {
        vec3 q = fibonacciPoint(ids[i], num);
        if (dot(q, nor) > 0.0) { // 外向きのみ
            testBump(ro, rd, sc, q, tmin, pos, col, occ);
        }
    }
}
    // 地面 y=0
    if (abs(rd.y) > EPS) {
        float h = (0.0 - ro.y) / rd.y;
        if (h > 0.0 && h + relEps(tmin)< tmin) {
            tmin = h;
            pos = ro + h * rd;
            nor = vec3(0.0, 1.0, 0.0);
            vec3 di  = sc - pos;
            float l  = max(length(di), EPS);
            occ = clamp(1.0 - dot(nor, di / l) / (l * l), 0.0, 1.0);
            col = vec3(1.0);
        }
    }
    if (tmin >= 1e19) { // 同じスケールのしきい値に
        col = vec3(0.6, 0.7, 0.9) * (0.65 + 0.35 * vv.y);
        fragColor = vec4(sqrt(col), 1.0);
        return;
    }

    pos = ro + tmin * rd;
    col *= occ;
    col = mix(col, vec3(1.0), 0.15 * (1.0 - exp(-0.002 * tmin * tmin)));


    // ガンマ補正
    col = sqrt(col);
    fragColor = vec4(col, 1.0);
}