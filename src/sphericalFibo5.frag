precision highp float;

uniform float u_time;
uniform vec2  u_resolution;
uniform float radius;          // 小球の半径係数（全体の密度感に影響）
uniform float num;             // 球面上の小球の総数
uniform float an;              // カメラの回転角
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
// [Level 1] 小球群（親の球の周り）との交差判定
//---------------------------------------------------------------
vec4 traceSmallSpheres(vec3 ro, vec3 rd, float rSmall) {
    // バウンディングスフィア判定（全体を包む球）
    float tBound = sphIntersect(ro, rd, vec3(0.0, 1.0, 0.0), 1.0 + rSmall);
    if(tBound < 0.0) return vec4(-1.0);

    // 衝突点から、探すべきフィボナッチグリッドを推定
    vec3 p = normalize((ro + tBound * rd) - vec3(0.0, 1.0, 0.0));

    // inverseSFのロジック
    float k = max(2.0, floor(log2(num * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z)) / log2(kPhi + 1.0)));
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));
    vec2  ka = 2.0 * F / num;
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));
    mat2  iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);
    vec2  c  = floor(iB * vec2(atan(p.y, p.x), p.z - 1.0 + 1.0 / num));

    float minT = 10000.0;
    float hitID = -1.0;
    
    // 近傍4点探索
    for (int s = 0; s < 4; s++) {
        vec2 uv = vec2(s & 1, s >> 1);
        float id = clamp(dot(F, uv + c), 0.0, num - 1.0);
        
        // 小球の中心位置
        vec3 center = getFibPos(id, num) + vec3(0.0, 1.0, 0.0);
        float t = sphIntersect(ro, rd, center, rSmall);
        
        if (t > 0.0 && t < minT) {
            minT = t;
            hitID = id;
        }
    }

    if(minT > 9999.0) return vec4(-1.0);
    return vec4(minT, hitID, 0.0, 0.0);
}

//---------------------------------------------------------------
// [Level 2] 小小球群（特定された小球の周り）との交差判定
// centerSmall: 親となる小球の中心座標
// radiusSmall: 親となる小球の半径
//---------------------------------------------------------------
vec4 traceMicroSpheres(vec3 ro, vec3 rd, vec3 centerSmall, float radiusSmall) {
    // マイクロ球の設定（密度は親と同じ num を使用）
    float numMicro = num; 
    
    // マイクロ球の半径
    // 小球の表面を埋めるため、小球の半径に対して比率で計算
    float rMicro = (2.0 / sqrt(numMicro)) * radiusSmall; 

    // 1. まず「親の小球」のバウンディング（少し大きめ）に当たるかチェック
    // マイクロ球が飛び出ている分、半径を +rMicro して判定
    float tBase = sphIntersect(ro, rd, centerSmall, radiusSmall + rMicro);
    if(tBase < 0.0) return vec4(-1.0);

    // 2. その衝突点における「小球表面上のローカル座標（方向）」を求める
    vec3 p = normalize((ro + tBase * rd) - centerSmall);

    // 3. ローカル座標を使ってフィボナッチグリッド推定
    float k = max(2.0, floor(log2(numMicro * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z)) / log2(kPhi + 1.0)));
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));
    vec2  ka = 2.0 * F / numMicro;
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));
    mat2  iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);
    vec2  c  = floor(iB * vec2(atan(p.y, p.x), p.z - 1.0 + 1.0 / numMicro));

    float minT = 10000.0;
    float hitID = -1.0;
    
    // 4. 近傍4点のマイクロ球との交差判定
    for (int s = 0; s < 4; s++) {
        vec2 uv = vec2(s & 1, s >> 1);
        float id = clamp(dot(F, uv + c), 0.0, numMicro - 1.0);
        
        // ★ここが重要★
        // マイクロ球の中心位置 = 親球の中心 + (フィボナッチ方向 * 親球の半径)
        // これにより、マイクロ球の中心が親球の表面上に位置し、半分埋まって半分飛び出る形になる
        vec3 dir = getFibPos(id, numMicro);
        vec3 microCenter = centerSmall + dir * radiusSmall;
        
        float t = sphIntersect(ro, rd, microCenter, rMicro);
        
        if (t > 0.0 && t < minT) {
            minT = t;
            hitID = id;
        }
    }

    if(minT > 9999.0) return vec4(-1.0);
    return vec4(minT, hitID, 0.0, 0.0);
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

    

    // ★ Step 1: 小球群（Level 1）との交差判定 ★
    // まず、どの「小球」のエリアにレイが飛んでいるか特定します
    vec4 res1 = traceSmallSpheres(ro, rd, sphereRad);
    
    // 小球のバウンディングエリアにヒットした場合
    if (res1.x > 0.0 && res1.x < tMin) {
        
        // ヒットした小球のIDと中心位置を取得
        float tSmallEst = res1.x; // これは小球そのものへの距離
        float idSmall = res1.y;
        vec3 centerSmall = getFibPos(idSmall, num) + vec3(0.0, 1.0, 0.0);
        
        // ★ Step 2: 小小球群（Level 2）との交差判定 ★
        // 特定した小球の表面に配置された「マイクロ球」との交差を詳しく調べる
        vec4 res2 = traceMicroSpheres(ro, rd, centerSmall, sphereRad);

        if (res2.x > 0.0 && res2.x < tMin) {
            // マイクロ球にヒット！
            tMin = res2.x;
            float idMicro = res2.y;
            
            // 法線の計算
            // 親と同じ密度(num)で配置されていると仮定して中心を再計算
            float numMicro = num; 
            vec3 dirMicro = getFibPos(idMicro, numMicro);
            
            // マイクロ球の中心
            // 親球の半径に「乗っかる」位置にある
            float rMicro = (2.0 / sqrt(numMicro)) * sphereRad;
            vec3 centerMicro = centerSmall + dirMicro * sphereRad;
            
            vec3 pos = ro + tMin * rd;
            nor = normalize(pos - centerMicro); // 法線はマイクロ球の中心から
            
            // 色: 親IDと子IDを混ぜて、カリフラワーのような複雑な色味を出す
            vec3 baseCol = 0.5 + 0.5 * sin(hash1(idSmall * 13.0) * 3.0 + 1.0 + vec3(0.0, 1.0, 1.0));
            // 少し明るさを変えて立体感を強調
            col = baseCol * (0.8 + 0.4 * sin(idMicro * 0.5));
            
        } else {
            // マイクロ球には当たらず、親の小球（の隙間）が見えている場合
            // ここを描画しないと隙間が透けてしまうので、親球の表面として描画
            if (tSmallEst < tMin) {
                tMin = tSmallEst;
                vec3 pos = ro + tMin * rd;
                nor = normalize(pos - centerSmall);
                
                // 隙間は少し暗くして、マイクロ球の出っ張りを強調する
                vec3 baseCol = 0.5 + 0.5 * sin(hash1(idSmall * 13.0) * 3.0 + 1.0 + vec3(0.0, 1.0, 1.0));
                col = baseCol * 0.5; 
            }
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