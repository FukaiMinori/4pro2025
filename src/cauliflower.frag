precision highp float;

uniform float uTime;
uniform vec2  u_resolution;
uniform float num;
uniform float an;
uniform float embedSmall;
uniform float embedSmall2;

out vec4 fragColor;

// Inigo Quilez - Golden spiral sphere mapping demo

const float kTau = 6.28318530718;            // 2π
const float kPhi = (1.0 + sqrt(5.0)) / 2.0;  // 黄金比 Φ
const float EPS  = 1e-6;

// 交差比較用の相対EPS
float relEps(float t) {
    return max(1e-6, 1e-5 * max(t, 1.0));
}

// 0~1 の乱数
float hash1(float n) { 
    return fract(sin(n) * 158.5453123);
}

// カリフラワー色
vec3 cauliflowerColor(float seed, float shade) {
    vec3 base = mix(
        vec3(0.95, 0.96, 0.90),       // 白
        vec3(0.8039, 0.9059, 0.7373), // 黄緑
        shade
    );

    float n = hash1(seed * 17.3);
    base += 0.05 * (n - 0.5);

    return clamp(base, 0.0, 1.0);
}

// 球との交差
bool intersectSphere(vec3 ro, vec3 rd, vec3 center, float radius, out float t) {
    vec3 oc = ro - center;
    float b = dot(rd, oc);
    float c = dot(oc, oc) - radius * radius;
    float D = b * b - c;
    if (D < 0.0) return false;
    float s  = sqrt(max(D, 0.0));
    float t0 = -b - s;
    float t1 = -b + s;
    t        = (t0 > 0.0) ? t0 : ((t1 > 0.0) ? t1 : -1.0);
    return t > 0.0;
}

// フィボナッチ関連
// 法線ベクトルから最近傍フィボナッチ点の index と距離
vec2 inverseSF(vec3 p) { 
    float k  = max(
        2.0, 
        floor(
            log2(num * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z))
          / log2(kPhi + 1.0)
        )
    );

    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));   // (F(k), F(k+1))

    vec2 ka = 2.0 * F / num;
    vec2 kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));

    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);

    vec2 c = floor(iB * vec2(
        atan(p.y, p.x),
        p.z - 1.0 + 1.0 / num
    ));

    float d = 8.0;
    float j = 0.0;

    for (int s = 0; s < 4; s++) {
        vec2  uv = vec2(s & 1, s >> 1);
        float id = clamp(dot(F, uv + c), 0.0, num - 1.0);

        float phi      = kTau * fract(id * kPhi);
        float cosTheta = 1.0 - (2.0 * id + 1.0) / num;
        float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

        vec3 q = vec3(
            cos(phi) * sinTheta,
            sin(phi) * sinTheta,
            cosTheta
        );
        float tmp = dot(q - p, q - p);
        if (tmp < d) {
            d = tmp;
            j = id;
        }
    }
    return vec2(j, sqrt(d));
}

// index → 単位球面上の点
vec3 fibonacciPoint(float id, float total) {
    float phi      = kTau * fract(id * kPhi);
    float cosTheta = 1.0 - (2.0 * id + 1.0) / total;
    float sinTheta = sqrt(max(0.0, 1.0 - cosTheta * cosTheta));
    return vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

// 近傍4候補の id を返す
void getNearest4Ids(vec3 p, out float ids[4]) {
    float k = max(
        2.0,
        floor(
            log2(num * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z))
          / log2(kPhi + 1.0)
        )
    );

    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));
    vec2 ka  = 2.0 * F / num;
    vec2 kb  = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));

    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);

    vec2 c = floor(iB * vec2(
        atan(p.y, p.x),
        p.z - 1.0 + 1.0 / num
    ));

    for (int s = 0; s < 4; s++) {
        vec2  uv = vec2(s & 1, s >> 1);
        float id = clamp(dot(F, uv + c), 0.0, num - 1.0);
        ids[s] = id;
    }
}

// カリフラワーの突起
// 大球外に突き出しているか
bool protrudesOut(vec3 center, float r, vec3 bigCenter, float bigR) {
    return (length(center - bigCenter) + r) > (bigR + 1e-5);
}

// 小小球ヒット共通処理
bool hitSubSphere(vec3 ro, vec3 rd, vec3 center, float radius,
    inout float tmin, inout vec3 pos, inout vec3 col, inout float occ
){
    float t;
    if (!intersectSphere(ro, rd, center, radius, t)) return false;
    if (t <= 0.0) return false;

    // 一番手前か？
    if (t >= tmin - relEps(tmin)) return false;

    vec3 hitPos = ro + t * rd;
    vec3 hitNor = normalize(hitPos - center);
    vec2 fi     = inverseSF(hitNor);

    float shade = clamp(0.6 + 0.4 * fi.y, 0.0, 1.0);
    vec3  c     = cauliflowerColor(fi.x + 10.0, shade);
    c *= 1.05;

    float o = 0.5 + 0.5 * hitNor.y;

    tmin = t;
    pos  = hitPos;
    col  = c;
    occ  = mix(1.0, o, 0.35);

    return true;
}

// 小球ヒット処理
bool hitParentSphere(vec3 ro, vec3 rd, vec3 center, float radius,
    inout float tmin, inout vec3 pos, inout vec3 col, inout float occ,
    out vec3 hitNor, out vec2 fiHit
){
    float tS;
    if (!intersectSphere(ro, rd, center, radius, tS)) return false;
    if (tS <= 0.0) return false;
    if (tS >= tmin - relEps(tmin)) return false; // ここだけで十分

    vec3 hitPos = ro + tS * rd;
    hitNor      = normalize(hitPos - center);
    fiHit       = inverseSF(hitNor);

    float shade   = clamp(0.4 + 0.6 * fiHit.y, 0.0, 1.0);
    vec3  iterCol = cauliflowerColor(fiHit.x, shade);
    iterCol      *= 0.90 + 0.20 * hitNor.y;
    float iterOcc = 0.2 + 0.2 * hitNor.y;

    tmin = tS;
    pos  = hitPos;
    col  = iterCol;
    occ  = mix(1.0, iterOcc, 0.35);

    return true;
}

// 一つのフィボナッチ点 q まわりの突起処理
void testBump(
    vec3 ro, vec3 rd, vec3 sc, vec3 q,
    inout float tmin, inout vec3 pos, inout vec3 col, inout float occ
){
    float bumpR0    = 0.3;
    float subRScale = 0.35;

    float rS      = bumpR0;
    float dCenter = 1.0 + rS * (1.0 - embedSmall);
    vec3  sphC    = sc + q * dCenter;
    float subR    = bumpR0 * subRScale;

    // 小小球（初期方向）
    vec3 q2_init   = normalize(q);
    vec3 subC_init = sphC + q2_init * (rS + subR * (1.0 - embedSmall2));
    hitSubSphere(ro, rd, subC_init, subR, tmin, pos, col, occ);

    // 小球
    vec3 hitNor;
    vec2 fiHit;
    if (hitParentSphere(ro, rd, sphC, rS, tmin, pos, col, occ, hitNor, fiHit))
    {
        // 小球ヒット後の精密化小小球
        vec3 q2_ref   = normalize(fibonacciPoint(fiHit.x, num));
        vec3 subC_ref = sphC + q2_ref * (rS + subR * (1.0 - embedSmall2));
        hitSubSphere(ro, rd, subC_ref, subR, tmin, pos, col, occ);
    }
}

// --------------------------------------------------------
// main
// --------------------------------------------------------

void main() {
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 p = (-u_resolution.xy + 2.0 * fragCoord) / u_resolution.y;

    // カメラ
    vec3 ro = vec3(2.5 * cos(an), 1.0, 2.5 * sin(an));
    vec3 ta = vec3(0.0, 1.0, 0.0);
    vec3 ww = normalize(ta - ro);
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0)));
    vec3 vv = normalize(cross(uu, ww));
    vec3 rd = normalize(p.x * uu + p.y * vv + 1.5 * ww);

    // 初期化
    vec3  sc   = vec3(0.0, 1.0, 0.0);   // 大球中心
    vec3  col  = vec3(0.0);
    float tmin = 1e20;
    float occ  = 1.0;
    vec3  pos  = vec3(0.0);
    vec3  nor  = vec3(0.0);
    vec2  fi   = vec2(0.0);

    // 大球ヒット
    float hBig;
    bool  hitBig = false;
    if (intersectSphere(ro, rd, sc, 1.0, hBig)) {
        if (hBig > 0.0 && hBig + relEps(tmin) < tmin) {
            tmin = hBig;
            pos  = ro + hBig * rd;
            nor  = normalize(pos - sc);
            occ  = 0.5 + 0.5 * nor.y;
            col  = vec3(0.92, 0.94, 0.88) * occ;
            fi   = inverseSF(nor);
            hitBig = true;
        }
    }

    // 背景側: 大球に当たらなかった場合の突起
    if (!hitBig) {
        float tBigOcc;
        bool  bigBlocks = intersectSphere(ro, rd, sc, 1.0, tBigOcc);

        vec3 pClosest = ro + rd * max(-dot(rd, ro - sc), 0.0);
        vec3 probeDir = normalize(pClosest - sc);

        float ids[4];
        getNearest4Ids(probeDir, ids);
        for (int i = 0; i < 4; i++) {
            vec3 q = fibonacciPoint(ids[i], num);

            float bumpR0    = 0.3;
            float subRScale = 0.35;
            float rS        = bumpR0;
            float dCenter   = 1.0 + rS * (1.0 - embedSmall);
            vec3 sphC       = sc + q * dCenter;
            float subR      = bumpR0 * subRScale;

            vec3 q2_init   = normalize(q);
            vec3 subC_init = sphC + q2_init * (rS + subR * (1.0 - embedSmall2));

            if (dot(q, probeDir) > 0.0 &&
                protrudesOut(subC_init, subR, sc, 1.0))
            {
                testBump(ro, rd, sc, q, tmin, pos, col, occ);
            }
        }
    }

    // 大球表面側: 大球に当たった場合の突起
    if (hitBig) {
        float ids[4];
        getNearest4Ids(nor, ids);
        for (int i = 0; i < 4; i++) {
            vec3 q = fibonacciPoint(ids[i], num);
            if (dot(q, nor) > 0.0) {
                testBump(ro, rd, sc, q, tmin, pos, col, occ);
            }
        }
    }

    // 地面 y = 0
    if (abs(rd.y) > EPS) {
        float h = (0.0 - ro.y) / rd.y;
        if (h > 0.0 && h + relEps(tmin) < tmin) {
            tmin = h;
            pos  = ro + h * rd;
            nor  = vec3(0.0, 1.0, 0.0);

            vec3  di = sc - pos;
            float l  = max(length(di), EPS);
            occ = clamp(1.0 - dot(nor, di / l) / (l * l), 0.0, 1.0);
            col = vec3(1.0);
        }
    }

    // どこにも当たらなかった場合の背景
    if (tmin >= 1e19) {
        fragColor = vec4(1.0);
        return;
    }
    pos = ro + tmin * rd;
    col *= occ;
    col = mix(col, vec3(1.0), 0.15 * (1.0 - exp(-0.002 * tmin * tmin)));

    // ガンマ補正
    col = sqrt(col);
    fragColor = vec4(col, 1.0);
}
