// golden_sphere_with_pearl.glsl
// 改良点: inverseSFで見つけた黄金螺旋点(fi)の位置を使い、その位置を中心に小球（パール）を生成します。
// 小球はメイン球の外側へ少しオフセットして描画され、既存の球の輪郭で切られないようにしています。

precision highp float;

uniform float u_time;
uniform vec2  u_resolution;
uniform float radius;     // （既存の用途のまま）
uniform float num;        // 球面上の点の総数
uniform float an;         // カメラ回転角
uniform float smallRadius; // 小球の半径
uniform float smallBias;   // 小球をメイン球表面からどれだけ外側にオフセットするか

out vec4 fragColor;

const float kTau = 6.28318530718;
const float kPhi = (1.0 + sqrt(5.0)) / 2.0;

// 0~1 の乱数
float hash1(float n) {
    return fract(sin(n) * 158.5453123);
}

// inverseSF は与えられた法線方向 p に対し、
// 最も近い黄金螺旋点の index とその点までの球面上の距離を返します。
vec2 inverseSF( vec3 p ) {
    float k  = max(
        2.0,
        floor(
            log2(num * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z))
            / log2(kPhi + 1.0)
        )
    );

    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));

    vec2  ka = 2.0 * F / num;
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));

    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x)
             / (ka.y * kb.x - ka.x * kb.y);

    vec2 c = floor(iB * vec2(
        atan(p.y, p.x),
        p.z - 1.0 + 1.0 / num
    ));

    float d = 8.0;
    float j = 0.0;

    for (int s = 0; s < 4; s++) {
        vec2  uv = vec2(s & 1, s >> 1);
        float id = clamp(dot(F, uv + c), 0.0, num - 1.0);

        float phi = kTau * fract(id * kPhi);
        float cosTheta = 1.0 - (2.0 * id + 1.0) / num;
        float sinTheta = sqrt(max(0.0, 1.0 - cosTheta * cosTheta));

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

void main() {
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 p = (-u_resolution.xy + 2.0 * fragCoord.xy) / u_resolution.y;

    // camera
    vec3 ro = vec3(2.5 * cos(an), 1.0, 2.5 * sin(an));
    vec3 ta = vec3(0.0, 1.0, 0.0);

    vec3 ww = normalize(ta - ro);
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0)));
    vec3 vv = normalize(cross(uu, ww));
    vec3 rd = normalize(p.x * uu + p.y * vv + 1.5 * ww);

    vec3 sc = vec3(0.0, 1.0, 0.0); // main sphere center

    vec3 col = vec3(1.0);
    float tmin = 10000.0;
    vec3 nor = vec3(0.0);
    float occ = 1.0;
    vec3 pos = vec3(0.0);

    // ground intersection
    float h = (0.0 - ro.y) / rd.y;
    if (h > 0.0) {
        tmin = h;
        nor = vec3(0.0, 1.0, 0.0);
        pos = ro + h * rd;
        vec3 di = sc - pos;
        float l = length(di);
        occ = 1.0 - dot(nor, di / l) / (l * l);
        col = vec3(1.0);
    }

    // main sphere intersection
    vec3 ce = ro - sc;
    float b = dot(rd, ce);
    float c = dot(ce, ce) - 1.0;
    h = b * b - c;

    bool hitMain = false;
    vec2 fi = vec2(0.0);
    vec3 mainHitPos = vec3(0.0);

    if (h > 0.0) {
        float t = -b - sqrt(h);
        if (t < tmin) {
            hitMain = true;
            tmin = t;
            mainHitPos = ro + t * rd;
            nor = normalize(mainHitPos - sc);
            occ = 0.5 + 0.5 * nor.y;

            // find nearest golden-spiral point info
            fi = inverseSF(nor);

            // original main-col setup (slightly toned down because small pearl may replace)
            col = 0.5 + 0.5 * sin(
                hash1(fi.x * 13.0) * 3.0 + 1.0 + vec3(0.0, 1.0, 1.0)
            );
            col *= smoothstep(0.02, 0.03, fi.y);
            col *= mix(1.0, 1.0 - smoothstep(0.12, 0.125, fi.y),
                        smoothstep(-0.9, 0.1, sin(u_time)));
            col *= 1.0 + 0.1 * sin(100.0 * radius * fi.y);
            col *= 1.2;
        }
    }

    // --- Create small sphere (pearl) centered at the golden-spiral point found by fi ---
    // Reconstruct the point q on unit sphere from the index fi.x
    float id = fi.x;
    float phi = kTau * fract(id * kPhi);
    float cosTheta = 1.0 - (2.0 * id + 1.0) / num;
    float sinTheta = sqrt(max(0.0, 1.0 - cosTheta * cosTheta));
    vec3 q = vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);

    // position of small sphere center: slightly pushed out along q to avoid clipping with main sphere
    float bias = smallBias; // e.g. 0.06
    vec3 pearlCenter = sc + q * (1.0 + bias);
    float rPearl = smallRadius; // e.g. 0.06

    // Ray-sphere intersection for small sphere
    vec3 ce_s = ro - pearlCenter;
    float b_s = dot(rd, ce_s);
    float c_s = dot(ce_s, ce_s) - rPearl * rPearl;
    float D_s = b_s * b_s - c_s;

    if (D_s > 0.0) {
        float t_s = -b_s - sqrt(D_s);
        if (t_s > 0.0 && t_s < tmin) {
            // small sphere is the closest hit: shade it clearly
            tmin = t_s;
            pos = ro + t_s * rd;
            nor = normalize(pos - pearlCenter);

            // lighting for pearl: soft diffuse + sharp specular to make it crisp
            vec3 lightDir = normalize(vec3(0.5, 1.0, 0.3));
            vec3 viewDir = normalize(ro - pos);
            float diff = max(0.0, dot(nor, lightDir));
            float rim = pow(1.0 - max(0.0, dot(viewDir, nor)), 3.0);
            float spec = pow(max(0.0, dot(reflect(-lightDir, nor), viewDir)), 64.0);

            // base color: slightly varied using the same hash as before for color harmony
            vec3 base = 0.6 + 0.4 * sin(hash1(fi.x * 7.0) * 6.28 + vec3(0.0, 2.0, 4.0));
            // pearls are typically desaturated and bright; mix with white for sheen
            vec3 pearlCol = mix(base, vec3(1.0), 0.6);

            col = pearlCol * (0.5 + 0.9 * diff) + vec3(1.0) * (0.9 * spec + 0.25 * rim);
            // enhance contrast
            col = pow(col, vec3(0.9));

            // small subtle outline to make it crisp against the main sphere
            float edge = smoothstep(rPearl * 0.98, rPearl * 1.02, length(pos - pearlCenter));
            col = mix(col, vec3(1.0), 1.0 - edge * 0.15);
        }
    }

    // if we hit main geometry (ground or main sphere) before a pearl hit, keep that shading
    if (tmin < 100.0) {
        if (tmin != 10000.0) {
            // pos already set for pearl; if not, compute for main hit
            if (tmin != h) {
                pos = ro + tmin * rd;
            }
            col *= occ;
            col = mix(col, vec3(1.0), 1.0 - exp(-0.003 * tmin * tmin));
        }
    }

    // ガンマ補正
    col = sqrt(max(col, vec3(0.0)));
    fragColor = vec4(col, 1.0);
}
