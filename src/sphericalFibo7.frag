precision highp float;

uniform float u_time;
uniform vec2  u_resolution;
uniform float radius;          // 小球の半径係数
uniform float num;             // 球面上の小球の総数
uniform float an;              // カメラの回転角
uniform float subNum;      // 小小球の総数（各小球の表面に打つ点の数）
uniform float subRadius;   // 小小球の半径
uniform float subOffset;   // 小小球の中心を小球表面からどれだけ外へ押し出すか


out vec4 fragColor;

//---------------------------------------------------------------
// 定数
//---------------------------------------------------------------
const float kTau = 6.28318530718;
const float kPhi = (1.0 + sqrt(5.0)) / 2.0;

//---------------------------------------------------------------
// 指定されたIDのフィボナッチ点の位置(単位球面上のベクトル)を返す
//---------------------------------------------------------------
vec3 getFibPos(float id, float n) {
    float phi = kTau * fract(id * kPhi);
    float cosTheta = 1.0 - (2.0 * id + 1.0) / n;
    float sinTheta = sqrt(clamp(1.0 - cosTheta * cosTheta, 0.0, 1.0));
    return vec3(
        cos(phi) * sinTheta,
        sin(phi) * sinTheta,
        cosTheta
    );
}

//---------------------------------------------------------------
// 単位球面上の座標 p に対して、最も近いフィボナッチ点までの情報を返す
// p: 表面上の位置(法線), n: 点の総数
// 戻り値: x=ID, y=距離
//---------------------------------------------------------------
vec2 getFibDist(vec3 p, float n) {
    // 層（k）を推定
    float k = max(2.0, floor(log2(n * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z)) / log2(kPhi + 1.0)));
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));
    
    vec2  ka = 2.0 * F / n;
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));
    mat2  iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);
    vec2  c  = floor(iB * vec2(atan(p.y, p.x), p.z - 1.0 + 1.0 / n));

    float d = 8.0;
    float j = 0.0;

    for (int s = 0; s < 4; s++) {
        vec2 uv = vec2(s & 1, s >> 1);
        float id = clamp(dot(F, uv + c), 0.0, n - 1.0);
        
        vec3 q = getFibPos(id, n); // ヘルパー関数を使用
        
        float tmp = dot(q - p, q - p);
        if (tmp < d) {
            d = tmp;
            j = id;
        }
    }
    return vec2(j, sqrt(d));
}

//---------------------------------------------------------------
// レイと球の交差判定関数
//---------------------------------------------------------------
float sphIntersect(vec3 ro, vec3 rd, vec3 ce, float ra) {
    vec3 oc = ro - ce;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - ra * ra;
    float h = b * b - c;
    if(h < 0.0) return -1.0;
    return -b - sqrt(h);
}

//---------------------------------------------------------------
// フィボナッチ球体群との交差判定を行うメイン関数
//---------------------------------------------------------------

vec4 traceFibSpheres(vec3 ro, vec3 rd, float rSmall) {
    float tBound = sphIntersect(ro, rd, vec3(0.0, 1.0, 0.0), 1.0 + rSmall + subRadius + subOffset);
    if(tBound < 0.0) return vec4(-1.0);

    vec3 p = normalize((ro + tBound * rd) - vec3(0.0, 1.0, 0.0));

    float k  = max(2.0, floor(log2(num * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z)) / log2(kPhi + 1.0)));
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));
    vec2  ka = 2.0 * F / num;
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));
    mat2  iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);
    vec2  c  = floor(iB * vec2(atan(p.y, p.x), p.z - 1.0 + 1.0 / num));

    float minT   = 10000.0;
    float hitID  = -1.0;
    float hitSub = 0.0;     // 0: 小球, 1: 小小球

    for (int s = 0; s < 4; s++) {
        vec2  uv     = vec2(s & 1, s >> 1);
        float id     = clamp(dot(F, uv + c), 0.0, num - 1.0);
        vec3  center = getFibPos(id, num) + vec3(0.0, 1.0, 0.0);

        // 小球の交差
        float tMain = sphIntersect(ro, rd, center, rSmall);
        if (tMain > 0.0 && tMain < minT) {
            minT  = tMain;
            hitID = id;
            hitSub = 0.0;
        }

        // 小小球の配置：視線が向かう局所方向から最近フィボナッチ点を推定
        // 近似方向（中心からレイの最近接点方向）
        vec3 probeDir = normalize((ro + max(tBound, 0.0) * rd) - center);
        if (subNum > 0.0 && subRadius > 0.0) {
            vec2 subInfo = getFibDist(probeDir, subNum);
            float subId  = subInfo.x;
            vec3 subDir  = getFibPos(subId, subNum);              // 単位球面上の方向
            vec3 subC    = center + subDir * (rSmall + subOffset + subRadius);

            float tSub = sphIntersect(ro, rd, subC, subRadius);
            if (tSub > 0.0 && tSub < minT) {
                minT   = tSub;
                hitID  = id;       // どの小球上の小小球かを保持
                hitSub = 1.0;      // 小小球ヒット
            }
        }
    }

    if(minT > 9999.0) return vec4(-1.0);
    // z: hitSub フラグ（0/1）
    return vec4(minT, hitID, hitSub, 0.0);
}

// ハッシュ関数
float hash1(float n) { 
    return fract(sin(n) * 158.5453123);
}

//---------------------------------------------------------------
// Main
//---------------------------------------------------------------
void main() {
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 p = (-u_resolution.xy + 2.0 * fragCoord.xy) / u_resolution.y;

    // カメラ設定
    vec3 ro = vec3(2.5 * cos(an), 1.0, 2.5 * sin(an));
    vec3 ta = vec3(0.0, 1.0, 0.0);
    vec3 ww = normalize(ta - ro);
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0)));
    vec3 vv = normalize(cross(uu, ww));
    vec3 rd = normalize(p.x * uu + p.y * vv + 1.5 * ww);

    vec3 col = vec3(1.0); // 背景色
    
    // 小球の半径計算
    float rScale = (radius > 0.0) ? radius : 1.0;
    float sphereRad = (2.0 / sqrt(num)) * rScale;

    // 地面との交差判定
    float tMin = 10000.0;
    vec3  nor = vec3(0.0);
    float hitID = -1.0;
    
    

    float hPlane = (0.0 - ro.y) / rd.y;
    if (hPlane > 0.0) {
        tMin = hPlane;
        nor = vec3(0.0, 1.0, 0.0);
        vec3 pos = ro + tMin * rd;
        vec3 di = vec3(0.0, 1.0, 0.0) - pos;
        float l = length(di);
        float occ = 1.0 - dot(nor, di / l) / (l * l);
        col = vec3(0.5) * occ;
    }

    // 球体群の判定
   // 球体群の判定
vec4 res = traceFibSpheres(ro, rd, sphereRad);

if (res.x > 0.0 && res.x < tMin) {
    tMin = res.x;
    hitID = res.y;

    vec3 pos = ro + tMin * rd;
    vec3 center = getFibPos(hitID, num) + vec3(0.0, 1.0, 0.0);

    // ベースカラー（小球ごとに色変化）
    vec3 baseCol = 0.5 + 0.5 * sin(hash1(hitID * 13.0) * 3.0 + 1.0 + vec3(0.0, 1.0, 1.0));

    if (res.z < 0.5) {
        // 小球ヒット
        nor = normalize(pos - center);
        col = baseCol;
    } else {
        // 小小球ヒット
        // もう一度小小球の中心を復元（trace と同じロジック）
        vec3 probeDir = normalize((ro + tMin * rd) - center);
        vec2 subInfo = getFibDist(probeDir, subNum);
        float subId  = subInfo.x;
        vec3 subDir  = getFibPos(subId, subNum);
        vec3 subC    = center + subDir * (sphereRad + subOffset + subRadius);

        nor = normalize(pos - subC);

        // 小小球の色：白ベース（出っ張りの白を継承）
        col = vec3(1.0);
    }
}


    // ライティング
    if(tMin < 100.0) {
        vec3 lig = normalize(vec3(1.0, 0.8, 0.6));
        float diff = clamp(dot(nor, lig), 0.0, 1.0);
        float amb = 0.5 + 0.5 * nor.y;
        col *= vec3(0.2, 0.3, 0.4) * amb + vec3(1.0, 0.9, 0.7) * diff;
        col = mix(col, vec3(1.0), 1.0 - exp(-0.003 * tMin * tMin));
    }

    col = sqrt(col);
    fragColor = vec4(col, 1.0);
}
