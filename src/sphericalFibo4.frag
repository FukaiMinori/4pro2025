precision highp float;

uniform float u_time;
uniform vec2  u_resolution;
uniform float radius; // 今回は未使用（内部で制御）
uniform float num;    // 500.0 ~ 1000.0 くらいが推奨
uniform float an;
out vec4 fragColor;

const float kTau = 6.28318530718;
const float kPhi = (1.0 + sqrt(5.0)) / 2.0;

//---------------------------------------------------------------
// ユーティリティ
//---------------------------------------------------------------
float hash1(float n) { return fract(sin(n) * 158.5453123); }

// スムーズな最小値（有機的な結合を作る）
float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    return mix(b, a, h) - k * h * (1.0 - h);
}

//---------------------------------------------------------------
// 黄金螺旋探索 (座標とIDを返す)
//---------------------------------------------------------------
vec4 inverseSF(vec3 p) {
    // 密度補正: 拡大縮小に対応するため num は固定で使用
    float localNum = num; 
    
    float k = max(2.0, floor(log2(localNum * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z)) / log2(kPhi + 1.0)));
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));
    vec2 ka = 2.0 * F / localNum;
    vec2 kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));
    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);
    vec2 c = floor(iB * vec2(atan(p.y, p.x), p.z - 1.0 + 1.0 / localNum));

    float d = 8.0;
    float bestID = 0.0;
    vec3  bestQ = vec3(0);

    for (int s = 0; s < 4; s++) {
        vec2 uv = vec2(s & 1, s >> 1);
        float id = clamp(dot(F, uv + c), 0.0, localNum - 1.0);
        float phi = kTau * fract(id * kPhi);
        float cosTheta = 1.0 - (2.0 * id + 1.0) / localNum;
        float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
        vec3 q = vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
        float tmp = dot(q - p, q - p);
        if (tmp < d) { d = tmp; bestID = id; bestQ = q; }
    }
    return vec4(bestQ, bestID);
}

//---------------------------------------------------------------
// 距離関数 (SDF) - フラクタル構造
//---------------------------------------------------------------
// 戻り値: x=距離, y=反復回数やIDを混ぜた色係数
vec2 map(vec3 p) {
    float d = length(p) - 1.0; // ベースの大球
    
    // パラメータ
    float scale = 1.0;         // 現在の空間スケール
    float orbit = 0.0;         // 色付け用の蓄積値
    
    vec3 currP = p;            // 変換されていく座標
    
    // 3回の反復（ループの中にループの構造）
    // 大球(Surface) -> 小球(Surface) -> 孫球(Surface)
    for(int i = 0; i < 3; i++) {
        // 1. 現在の座標空間における球表面での「最寄り点(fi)」を探す
        // 原点からの方向ベクトルを使用
        vec3 dir = normalize(currP);
        vec4 res = inverseSF(dir); // xyz: 中心, w: ID
        vec3 center = res.xyz;

        // 2. 空間の折り畳み (Folding)
        // 座標系を「見つけた小球の中心」へ移動させる
        // 半径1.0の表面にある点へ移動
        currP = currP - center * 1.0;
        
        // 3. 次の階層のために空間を拡大 (Zoom in)
        // 小球のサイズ比率 (密度 num に依存して隙間なく埋める係数を設定)
        // 3.0 / sqrt(num) は「小球が隣と接するくらい」の経験則
        float sizeRatio = 2.8 / sqrt(num);
        
        // 拡大
        currP /= sizeRatio;
        scale *= sizeRatio;
        
        // 4. 距離の合成
        // 現在の空間での「単位球」との距離を計算し、全体スケールで割って実距離に戻す
        float distSphere = (length(currP) - 1.0) * scale;
        
        // 有機的に結合 (Smooth Min)
        // 階層が深くなるほど結合を鋭くする (* scale)
        d = smin(d, distSphere, 0.15 * scale);
        
        // 色情報の蓄積
        orbit += hash1(res.w + float(i)*13.52) * 0.5;
    }

    return vec2(d, orbit);
}

//---------------------------------------------------------------
// 法線計算
//---------------------------------------------------------------
vec3 calcNormal(vec3 p) {
    vec2 e = vec2(0.0005, 0.0); // 精度向上のため少し小さく
    return normalize(vec3(
        map(p + e.xyy).x - map(p - e.xyy).x,
        map(p + e.yxy).x - map(p - e.yxy).x,
        map(p + e.yyx).x - map(p - e.yyx).x
    ));
}

//---------------------------------------------------------------
// Shadow (簡易ソフトシャドウ)
//---------------------------------------------------------------
float softShadow(vec3 ro, vec3 rd, float k) {
    float res = 1.0;
    float t = 0.02;
    for(int i=0; i<12; i++) {
        float h = map(ro + rd*t).x;
        res = min(res, k*h/t);
        t += clamp(h, 0.02, 0.2);
        if(res < 0.001 || t > 2.0) break;
    }
    return clamp(res, 0.0, 1.0);
}

//---------------------------------------------------------------
// main
//---------------------------------------------------------------
void main() {
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 p = (-u_resolution.xy + 2.0 * fragCoord.xy) / u_resolution.y;

    // カメラ回転
    vec3 ro = vec3(2.8 * cos(an), 1.0, 2.8 * sin(an)); // 少し引いて全体を見る
    vec3 ta = vec3(0.0, 0.0, 0.0);
    vec3 ww = normalize(ta - ro);
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0)));
    vec3 vv = normalize(cross(uu, ww));
    vec3 rd = normalize(p.x * uu + p.y * vv + 1.5 * ww);

    vec3 col = vec3(0.02, 0.03, 0.05); // 暗い背景
    
    // レイマーチング
    float t = 0.0;
    float tmax = 10.0;
    vec2 res = vec2(0.0);
    bool hit = false;
    
    for(int i=0; i<90; i++) {
        vec3 pos = ro + t * rd;
        res = map(pos); // x:距離, y:色ID
        if(res.x < 0.001) {
            hit = true;
            break;
        }
        t += res.x * 0.6; // 安全のためステップ幅を少し抑える
        if(t > tmax) break;
    }

    if(hit) {
        vec3 pos = ro + t * rd;
        vec3 nor = calcNormal(pos);
        
        // 色の生成 (蓄積されたorbit値を使用)
        float idVal = res.y;
        vec3 baseCol = 0.5 + 0.5 * sin(vec3(0.0, 0.3, 0.6) + idVal * 3.0);
        // カリフラワーっぽく少し白/クリーム色を混ぜる
        baseCol = mix(baseCol, vec3(0.9, 0.85, 0.7), 0.5);

        // ライティング
        vec3  lig = normalize(vec3(1.0, 0.8, 0.6));
        float dif = clamp(dot(nor, lig), 0.0, 1.0);
        float amb = 0.5 + 0.5 * nor.y;
        float sha = softShadow(pos + nor*0.01, lig, 8.0); // 自身の影
        
        // SSS (Subsurface Scattering) 風の透け感
        // 法線と視線、ライトの関係で「薄い部分」を明るく見せるハック
        float sss = pow(clamp(1.0 + dot(nor, rd), 0.0, 1.0), 2.0) * 0.5;

        vec3 lin = vec3(0.0);
        lin += 1.5 * dif * vec3(1.0, 0.95, 0.8) * sha;
        lin += 0.4 * amb * vec3(0.3, 0.4, 0.6); // 青っぽい環境光
        lin += 1.0 * sss * vec3(1.0, 0.4, 0.2); // 内部散乱光（赤み）

        col = baseCol * lin;
        
        // フォグ (奥を暗く)
        col = mix(col, vec3(0.02, 0.03, 0.05), 1.0 - exp(-0.1 * t * t));
    }

    // ガンマ補正
    col = pow(col, vec3(1.0/2.2));
    fragColor = vec4(col, 1.0);
}