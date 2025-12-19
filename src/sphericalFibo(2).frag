#version 300 es
precision highp float;

uniform float u_time;
uniform vec2  u_resolution;
out vec4 fragColor;

// ---------------------------------------------------------------
// Inigo Quilez - Golden spiral sphere mapping demo (元コードを拡張)
// ---------------------------------------------------------------

vec2 inverseSF( vec3 p ) { 
    const float kTau = 6.28318530718;
    const float kPhi = (1.0 + sqrt(5.0)) / 2.0;
    const float kNum = 130.0;

    float k  = max(2.0, floor(log2(kNum * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z)) / log2(kPhi + 1.0)));
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));

    vec2  ka = 2.0 * F / kNum;
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));

    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);
    vec2 c = floor(iB * vec2(atan(p.y, p.x), p.z - 1.0 + 1.0 / kNum));

    float d = 8.0;
    float j = 0.0;
    for (int s = 0; s < 4; s++) {
        vec2  uv = vec2(s & 1, s >> 1);
        float id = clamp(dot(F, uv + c), 0.0, kNum - 1.0);

        float phi = kTau * fract(id * kPhi);
        float cosTheta = 1.0 - (2.0 * id + 1.0) / kNum;
        float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

        vec3 q = vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
        float tmp = dot(q - p, q - p);
        if (tmp < d) {
            d = tmp;
            j = id;
        }
    }
    return vec2(j, sqrt(d));
}

float hash1(float n) { return fract(sin(n) * 158.5453123); }

// --- ヘルパ：フィボナッチ（黄金螺旋）上の i 番目点（単位球上）を返す
vec3 fiboPoint(float id) {
    const float kTau = 6.28318530718;
    const float kPhi = (1.0 + sqrt(5.0)) / 2.0;
    const float kNum = 130.0;
    float phi = kTau * fract(id * kPhi);
    float cosTheta = 1.0 - (2.0 * id + 1.0) / kNum;
    float sinTheta = sqrt(max(0.0, 1.0 - cosTheta * cosTheta));
    return vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

// --- レイと球の解析的交差（返り値：最小正の t、なければ -1）
// center 'c', 半径 'r'
float intersectSphereAnalytic(vec3 ro, vec3 rd, vec3 c, float r) {
    vec3 oc = ro - c;
    float b = dot(rd, oc);
    float c0 = dot(oc, oc) - r * r;
    float disc = b * b - c0;
    if (disc <= 0.0) return -1.0;
    float t = -b - sqrt(disc);
    if (t > 0.0001) return t;
    t = -b + sqrt(disc);
    if (t > 0.0001) return t;
    return -1.0;
}

void main() {
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 p = (-u_resolution.xy + 2.0 * fragCoord.xy) / u_resolution.y;

    // カメラ
    float an = 0.0; // 回転いれるなら 0.5 * u_time など
    vec3 ro = vec3(2.5 * cos(an), 1.0, 2.5 * sin(an));
    vec3 ta = vec3(0.0, 1.0, 0.0);
    vec3 ww = normalize(ta - ro);
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0)));
    vec3 vv = normalize(cross(uu, ww));
    vec3 rd = normalize(p.x * uu + p.y * vv + 1.5 * ww);

    // シーン要素
    vec3 sc = vec3(0.0, 1.0, 0.0); // 大球の中心
    float tmin = 1e6;
    vec3  nor = vec3(0.0);
    vec3  pos = vec3(0.0);
    float occ = 1.0;

    vec3 col = vec3(1.0);

    // 地面 y=0 と交差
    float h = (0.0 - ro.y) / rd.y;
    if (h > 0.0) {
        tmin = h;
        nor = vec3(0.0, 1.0, 0.0);
        pos = ro + h * rd;
        vec3 di = sc - pos;
        float l = length(di);
        occ = 1.0 - dot(nor, di / l) * 1.0 / (l * l);
        col = vec3(1.0);
    }

    // -------------------------
    // 大球（解析的交差、半径=1）
    // -------------------------
    {
        vec3 ce = ro - sc;
        float b = dot(rd, ce);
        float c = dot(ce, ce) - 1.0;
        float disc = b * b - c;
        if (disc > 0.0) {
            float tSphere = -b - sqrt(disc);
            if (tSphere > 0.0001 && tSphere < tmin) {
                tmin = tSphere;
                pos = ro + tmin * rd;
                nor = normalize(pos - sc);
                occ = 0.5 + 0.5 * nor.y;

                vec2 fi = inverseSF(nor); // fi.x = index, fi.y = 距離
                col = 0.5 + 0.5 * sin(hash1(fi.x * 13.0) * 3.0 + 1.0 + vec3(0.0, 1.0, 1.0));
                col *= smoothstep(0.02, 0.03, fi.y);
                col *= mix(1.0, 1.0 - smoothstep(0.12, 0.125, fi.y), smoothstep(-0.1, 0.1, sin(u_time)));
                col *= 1.0 + 0.1 * sin(250.0 * fi.y);
                col *= 1.5;
            }
        }
    }

    
    const float rSmall = 0.2; // 小球半径
    const int K = 130;
    for (int i = 0; i < K; i++) {
        float id = float(i);
        vec3 q = fiboPoint(id); // 単位球上の点
        vec3 center = sc + q * (1.0 + rSmall); // 大球表面に接するように配置

        float ti = intersectSphereAnalytic(ro, rd, center, rSmall);
        if (ti > 0.0 && ti < tmin) {
            tmin = ti;
            pos = ro + tmin * rd;
            nor = normalize(pos - center);

            float hue = hash1(id * 13.0);
            col = 0.3 + 0.5 * sin(hue * 3.0 + 1.0 + vec3(0.0, 1.0, 1.0));
            float spec = pow(max(0.0, dot(normalize(vec3(0.4,0.8,0.2)), nor)), 16.0);
            col += vec3(spec * 0.6);
            occ = 0.5 + 0.5 * nor.y;
        }
    }

    // 背景・減衰処理
    if (tmin < 1e5) {
        pos = ro + tmin * rd;
        col *= occ;
        col = mix(col, vec3(1.0), 1.0 - exp(-0.003 * tmin * tmin));
    }

    // ガンマ補正と出力
    col = sqrt(max(vec3(0.0), col));
    fragColor = vec4(col, 1.0);
}
