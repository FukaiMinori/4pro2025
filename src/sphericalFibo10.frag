precision highp float;                 // 浮動小数点精度を高く設定（演算誤差を低減）

uniform float uTime;                  // 経過時間
uniform vec2  u_resolution;            // 画面解像度（幅, 高さ）
uniform float radius;                  // 半径（最終的に色のゆらぎで使用）
uniform float num;                     // 球面上の点の総数
uniform float an;                      // カメラの回転角
out vec4 fragColor;                    // 出力するピクセルの最終 RGBA 値

//---------------------------------------------------------------
// Inigo Quilez - Golden spiral sphere mapping demo
//---------------------------------------------------------------

// p: 球面上の法線ベクトル（単位ベクトル）
// 戻り値: 
//   x = 最も近い黄金螺旋点のインデックス
//   y = その点との距離
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

// 黄金螺旋点（インデックス→単位球座標）
vec3 fibonacciPoint(float id, float total) {
    const float kTau = 6.28318530718;
    const float kPhi = (1.0 + sqrt(5.0)) / 2.0;
    float phi = kTau * fract(id * kPhi);
    float cosTheta = 1.0 - (2.0 * id + 1.0) / total;
    float sinTheta = sqrt(max(0.0, 1.0 - cosTheta * cosTheta));
    return vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

vec2 hash22(vec2 p){
    float n = dot(p, vec2(127.1, 311.7));
    float m = dot(p, vec2(269.5, 183.3));
    return fract(sin(vec2(n, m)) * 43758.5453);
}

// 球面での微小ジッター：特徴点 q を接線方向に小回転
vec3 jitterOnSphere(vec3 q, float seed, float amp) {
    // 任意の接線基底（qと直交）
    vec3 t1 = normalize(abs(q.z) < 0.99 ? cross(q, vec3(0,0,1)) : cross(q, vec3(0,1,0)));
    vec3 t2 = normalize(cross(q, t1));
    // 角度と向きをハッシュで決定
    float a = 6.2831853 * fract(sin(seed * 13.37) * 157.31); // [0,2π)
    float r = amp * (fract(sin(seed * 91.7) * 43758.2) - 0.5); // [-amp/2, amp/2]
    vec3 dir = cos(a) * t1 + sin(a) * t2;
    // 小回転の一次近似（弦距離ベースならこの擾動で十分）
    vec3 qj = normalize(q + r * dir);
    return qj;
}

// 球面 Worley の F1/F2 を求める（黄金螺旋近傍 12 サンプル）
struct WorleyRes { float f1; float f2; float id1; float id2; };

WorleyRes worleySphereNearest2(vec3 p, float total, float jitterAmp) {
    const float kTau = 6.28318530718;
    const float kPhi = (1.0 + sqrt(5.0)) / 2.0;

    // 近傍層の推定（inverseSF と同じ）
    float k  = max(2.0, floor(
        log2(total * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z))
        / log2(kPhi + 1.0)
    ));

    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));

    vec2 ka = 2.0 * F / total;
    vec2 kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0));

    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x) / (ka.y * kb.x - ka.x * kb.y);

    vec2 c = floor(iB * vec2(atan(p.y, p.x), p.z - 1.0 + 1.0 / total));

    // サンプルオフセット（3×4=12点）
    vec2 OFFS[12] = vec2[12](
        vec2(-1.0,-1.0), vec2(0.0,-1.0), vec2(1.0,-1.0),
        vec2(-1.0, 0.0), vec2(0.0, 0.0), vec2(1.0, 0.0),
        vec2(-1.0, 1.0), vec2(0.0, 1.0), vec2(1.0, 1.0),
        vec2(2.0, 0.0),  vec2(0.0, 2.0), vec2(-2.0, 0.0)
    );

    float best1 = 1e9, best2 = 1e9;
    float id1 = -1.0, id2 = -1.0;

    for (int s = 0; s < 12; s++) {
        vec2 uv = OFFS[s] + c;
        float id = clamp(dot(F, uv), 0.0, total - 1.0);

        // 特徴点の位置（ジッターあり）
        vec3 q = fibonacciPoint(id, total);
        // ジッター量は弦距離スケールに合わせて小さめに（例：0.01〜0.05）
        q = jitterOnSphere(q, id + total, jitterAmp);

        // 弦距離（= |q - p|）
        float d = length(q - p);

        if (d < best1) {
            best2 = best1; id2 = id2;
            best1 = d;     id2 = id1;
            id1   = id;
        } else if (d < best2) {
            best2 = d;     id2 = id;
        }
    }

    WorleyRes wr;
    wr.f1 = best1; wr.f2 = best2;
    wr.id1 = id1;  wr.id2 = id2;
    return wr;
}

// 0~1 の乱数を返す簡易 hash 関数
float hash1(float n) { 
    return fract(sin(n) * 158.5453123);
}

//---------------------------------------------------------------
// main: ピクセルごとの描画処理
//---------------------------------------------------------------
void main() {

    // ピクセル座標を正規化（中心を (0,0) に）
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 p = (-u_resolution.xy + 2.0 * fragCoord.xy) / u_resolution.y;

    // カメラ設定
    //カメラ位置。anにより、カメラが円周上を回る。高さは1.0
    vec3 ro = vec3(2.5 * cos(an), 1.0, 2.5 * sin(an));
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


    //ro = r0 = レイの始点(カメラ位置)
    //rd = レイの方向
    //sc = 球の中心C(vec3(0.0,1.0,0.0))
    //t(一時変数hを利用している) = レイパラメータ(距離)
    //tmin = 現在見つかっている最小の交差距離(小さいほど手前で交差している)
    //pos = 交差点位置 = ro + t*rd
    //nor = 交差点の法線(球ならnormalize(pos - sc)、平面なら固定)

    // 地面 y=0 との交差判定
    
    float h = (0.0 - ro.y) / rd.y;
    //tを解く式(平面y=0と交差するt)
    //交差点位置でのレイのy成分：ro.y + trd.y = 0→この式をtについて解くと、hの式になる。
    if (h > 0.0) { //レイが前方(カメラから見て正方向)に交差しているかの確認
        tmin = h; //最初の交差距離
        nor = vec3(0.0, 1.0, 0.0); //地面の法線(上向き)
        pos = ro + h * rd; //交点座標の計算
        vec3 di = sc - pos; //地面上の交点から球中心までのベクトル(球との距離や影処理に使用)
        float l = length(di); //diの長さ
        occ = 1.0 - dot(nor, di / l) / (l * l);    // 影の処理
        col = vec3(1.0);
    }

    // 球との交差判定
    vec3  ce = ro - sc;
    float b = dot(rd, ce);
    float c = dot(ce, ce) - 1.0;   
    h = b * b - c;                 //判別式

    if (h > 0.0) {
        h = -b - sqrt(h);          //t1→レイが球と最初に当たる交差距離
        if (h < tmin) {
            tmin = h;
            nor = normalize(ro + h * rd - sc); // 法線
            occ = 0.5 + 0.5 * nor.y;           // yが大きい(上向きほど明るい
        }

        //球面の法線をinverseSFに渡し、その方向上に最も近い螺旋状の点のインデックスとその点までの距離を返す
        vec2 fi = inverseSF(nor);

        // Worley（F1/F2）
        WorleyRes wr = worleySphereNearest2(nor, num, 0.03); // ジッター振幅は調整可

        // 代表的なセルラーノイズ指標
        float F1 = wr.f1;            // 最近傍距離
        float F2 = wr.f2;            // 第2近傍距離
        float cell = F2 - F1;        // Worley のよく使うコンビネーション
        float ridges = 1.0 - smoothstep(0.0, 0.09, F1); // ボーダー強調例

        // 点ごとのカラーランダム（インデックス依存）
        /*
        float rnd = hash1(wr.id1 * 13.0);
        vec3 base = vec3(0.35, 0.6, 1.0);
        vec3 tint = mix(vec3(1.0), vec3(0.9, 0.95, 1.2), rnd);
        */

        // 色決定（例）
        /*
        col = base * tint;
        col *= 0.6 + 0.4 * smoothstep(0.02, 0.14, cell); // セル内部の明暗
        col *= 0.7 + 0.3 * ridges;                       // 境界リムを少し足す
        col *= 1.5;
        */

        //vec2 fi = inverseSF(nor);
        float phase = uTime + fi.x * 0.1;
        float twinkle = 0.5 + 0.5 * sin(phase);
        float glow = exp(-40.0 * fi.y * fi.y);
        
        float rnd = hash1(fi.x * 17.0);
        vec3 starColor = mix(vec3(0.8,0.9,1.0), vec3(1.0,0.8,0.6), rnd);
        
        col = starColor * (0.5 + glow * twinkle);
    }
    
    //---------------------------
    // 背景フェード
    //---------------------------
    if (tmin < 100.0) {
        pos = ro + tmin * rd;
        col *= occ;
        col = mix(col, vec3(1.0), 1.0 - exp(-0.003 * tmin * tmin));
    }

    //---------------------------
    // ガンマ補正
    //---------------------------
    col = sqrt(col);
    fragColor = vec4(col, 1.0);
}
