precision highp float;

uniform float u_time;
uniform vec2  u_resolution;
uniform float radius;    // 小球の半径係数 (元からの)
uniform float num;       // 球面上の小球の総数
uniform float an;        // カメラ回転角
uniform float subNum;    // 各小球に付く小小球の数 (例: 30.0)
uniform float tinyScale; // 小小球のサイズ係数 (例: 0.25)
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
// レイと球の交差判定関数 (正・負のtを返さない)
// ro: ray origin, rd: ray dir, ce: center, ra: radius
// 戻り値: 最小正解 t (ヒットがなければ -1.0)
//---------------------------------------------------------------
float sphIntersect(vec3 ro, vec3 rd, vec3 ce, float ra) {
    vec3 oc = ro - ce;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - ra * ra;
    float h = b * b - c;
    if (h < 0.0) return -1.0;
    float t = -b - sqrt(h);
    if (t > 0.0001) return t; // 近接交差を除去
    t = -b + sqrt(h);
    if (t > 0.0001) return t;
    return -1.0;
}

//---------------------------------------------------------------
// フィボナッチ球体群との交差判定（小球＋各小球の小小球も調べる）
// 戻り値: vec4(minT, smallID, subID, typeFlag)
//   minT: 距離
//   smallID: ヒットした「小球」のID（-1: 無し）
//   subID: ヒットした小小球のID（-1: 小球本体にヒット、>=0: 小小球のインデックス）
//   typeFlag: 0 = none, 1 = small sphere, 2 = tiny sphere
//---------------------------------------------------------------
vec4 traceFibSpheres(vec3 ro, vec3 rd, float rSmall, float subN, float tinyS) {
    // バウンディングスフィア判定（全体）
    float tBound = sphIntersect(ro, rd, vec3(0.0, 1.0, 0.0), 1.0 + rSmall + tinyS);
    if (tBound < 0.0) return vec4(-1.0, -1.0, -1.0, 0.0);

    // グリッド推定（既存ロジックを簡潔に再利用）
    vec3 p = normalize((ro + tBound * rd) - vec3(0.0, 1.0, 0.0));

    float k = max(2.0, floor(log2(num * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z)) / log2(kPhi + 1.0)));
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));
    vec2  ka = 2.0 * F / num;
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));
    mat2  iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);
    vec2  c  = floor(iB * vec2(atan(p.y, p.x), p.z - 1.0 + 1.0 / num));

    float minT = 1e5;
    float hitSmallID = -1.0;
    float hitSubID = -1.0;
    float hitType = 0.0;

    // 近傍4点探索（各候補小球に対して小球本体とその小小球群をチェック）
    for (int s = 0; s < 4; s++) {
        vec2 uv = vec2(float(s & 1), float(s >> 1));
        float id = clamp(dot(F, uv + c), 0.0, num - 1.0);

        vec3 center = getFibPos(id, num) + vec3(0.0, 1.0, 0.0);

        // 1) 小球本体との交差
        float tSmall = sphIntersect(ro, rd, center, rSmall);
        if (tSmall > 0.0 && tSmall < minT) {
            minT = tSmall;
            hitSmallID = id;
            hitSubID = -1.0; // 小球本体
            hitType = 1.0;
        }

        // 2) その小球に付く各小小球（フィボナッチ点ごと）
        // 注意: subN が float のため int に変換
        int nSub = int(clamp(subN, 0.0, 256.0)); // 上限で保護
        for (int si = 0; si < nSub; si++) {
            float fid = float(si);
            vec3 local = getFibPos(fid, subN); // 単位球上の位置（ローカル）
            // 小小球の中心は小球の表面上 (外向き)
            vec3 tinyCenter = center + local * rSmall;
            float tinyRad = rSmall * tinyS;
            float tTiny = sphIntersect(ro, rd, tinyCenter, tinyRad);
            if (tTiny > 0.0 && tTiny < minT) {
                minT = tTiny;
                hitSmallID = id;
                hitSubID = fid; // 小小球のインデックス
                hitType = 2.0;
            }
        }
    }

    if (minT > 90000.0) return vec4(-1.0, -1.0, -1.0, 0.0);
    return vec4(minT, hitSmallID, hitSubID, hitType);
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

    vec3 col = vec3(1.0); // 背景色（明るめ）
    
    // 小球（attached）の半径計算
    float rScale = (radius > 0.0) ? radius : 1.0;
    float sphereRad = (2.0 / sqrt(num)) * rScale;

    // 地面との交差判定（簡易）
    float tMin = 1e5;
    vec3 nor = vec3(0.0);
    float hitSmallID = -1.0; // ヒットした小球ID
    float hitSubID = -1.0;   // ヒットした小小球ID
    float hitType = 0.0;     // 1=small,2=tiny

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

    // 球体群（小球 + 小小球）との判定
    float tinyScaleClamped = clamp(tinyScale, 0.01, 1.0);
    vec4 res = traceFibSpheres(ro, rd, sphereRad, subNum, tinyScaleClamped);

    if (res.x > 0.0 && res.x < tMin) {
        tMin = res.x;
        hitSmallID = res.y;
        hitSubID = res.z;
        hitType = res.w;

        // 衝突点
        vec3 pos = ro + tMin * rd;

        if (hitType == 1.0) {
            // 小球本体にヒット
            vec3 center = getFibPos(hitSmallID, num) + vec3(0.0, 1.0, 0.0);
            nor = normalize(pos - center);
            vec3 baseCol = 0.5 + 0.5 * sin(hash1(hitSmallID * 13.0) * 3.0 + 1.0 + vec3(0.0, 1.0, 1.0));
            // その小球上のフィボナッチ点を使って「ドット表示」も可能（既存の表現）
            float subPoints = subNum;
            vec2 subFib = vec2(0.0);
            if (subPoints > 0.5) {
                // nor は単位ベクトルなので、そのまま getFibPos inverse を使えるわけではないが
                // ここでは既存の getFibDist を使えないため、代わりに小さな光沢表現は省略
            }
            col = baseCol;
        } else if (hitType == 2.0) {
            // 小小球（tiny）にヒット
            // 小小球の中心と法線を再計算
            float smallId = hitSmallID;
            float subId = hitSubID;
            vec3 smallCenter = getFibPos(smallId, num) + vec3(0.0, 1.0, 0.0);
            vec3 local = getFibPos(subId, subNum); // 単位方向
            vec3 tinyCenter = smallCenter + local * sphereRad;
            vec3 tinyNormal = normalize(pos - tinyCenter);
            nor = tinyNormal;

            // 小小球は白くする（若干ハイライト）
            col = vec3(1.0);
        }
    }

    // ライティング（影や反射は簡易）
    if (tMin < 1e4) {
        vec3 lig = normalize(vec3(1.0, 0.8, 0.6));
        float diff = clamp(dot(nor, lig), 0.0, 1.0);
        float amb = 0.5 + 0.5 * nor.y;
        // 小小球は白（col already set）、小球は base色（col set）
        col *= vec3(0.2, 0.3, 0.4) * amb + vec3(1.0, 0.9, 0.7) * diff;
        // 距離減衰でフェード（シーンの深度感）
        col = mix(col, vec3(1.0), 1.0 - exp(-0.003 * tMin * tMin));
    }

    col = sqrt(max(col, vec3(0.0))); // ガンマっぽく
    fragColor = vec4(col, 1.0);
}
