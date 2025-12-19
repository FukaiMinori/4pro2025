precision highp float;                 // 浮動小数点精度を高く設定（演算誤差を低減）

uniform float u_time;                  // 経過時間
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

    //---------------------------
    // 層（k）を推定
    // 高さ p.z から黄金螺旋ライン上のどの層にいるかを対数式で近似
    //---------------------------
    float k  = max(
        2.0, 
        floor(
            log2(num * kTau * 0.5 * sqrt(5.0) * (1.0 - p.z * p.z))
            / log2(kPhi + 1.0)
        )
    );

    //---------------------------
    // フィボナッチ近似を使って F(k), F(k+1) を推定
    //---------------------------
    float Fk = pow(kPhi, k) / sqrt(5.0);
    vec2  F  = vec2(round(Fk), round(Fk * kPhi));   // (F(k), F(k+1))

    //---------------------------
    // 黄金螺旋グリッドの係数を生成
    //---------------------------
    vec2  ka = 2.0 * F / num;                       // 緯度方向の係数
    vec2  kb = kTau * (fract((F + 1.0) * kPhi) - (kPhi - 1.0)); // 経度方向の係数

    //---------------------------
    // 逆行列により点 p に最も近い格子番号 c を推定
    //---------------------------
    mat2 iB = mat2(ka.y, -ka.x, kb.y, -kb.x)
             / (ka.y * kb.x - ka.x * kb.y);

    vec2 c = floor(iB * vec2(
        atan(p.y, p.x),                  // 経度 φ
        p.z - 1.0 + 1.0 / num            // 緯度方向補正
    ));

    //---------------------------
    // 近傍の 4 点をチェックして最短距離の点を探す
    //---------------------------
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

//---------------------------------------------------------------
// 0~1 の乱数を返す簡易 hash 関数
//---------------------------------------------------------------
float hash1(float n) { 
    return fract(sin(n) * 158.5453123);
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

//---------------------------------------------------------------
// main: ピクセルごとの描画処理
//---------------------------------------------------------------
void main() {

    //---------------------------
    // ピクセル座標を正規化（中心を (0,0) に）
    //---------------------------
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 p = (-u_resolution.xy + 2.0 * fragCoord.xy) / u_resolution.y;

    //---------------------------
    // カメラ設定
    //---------------------------
    //カメラ位置。anにより、カメラが円周上を回る。高さは1.0
    vec3 ro = vec3(2.5 * cos(an), 1.0, 2.5 * sin(an));
    vec3 ta = vec3(0.0, 1.0, 0.0); // 注視点。カメラはここに向かっている

    vec3 ww = normalize(ta - ro); //カメラの前方向(正規化したta-ro)
    vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0))); //右方向ベクトル(wwとy軸との外積)
    vec3 vv = normalize(cross(uu, ww));  //上方向ベクトル(uuとwwの外積)
    //スクリーン上の座標p.xとp.yをカメラの右・上ベクトルに乗せ、前方向に1.5*wwを加えることで視野を作る
    //1.5は視野に相当する値
    vec3 rd = normalize(p.x * uu + p.y * vv + 1.5 * ww); 

    //---------------------------
    // 初期化
    //---------------------------
    vec3 sc = vec3(0.0, 1.0, 0.0); // 球の中心
    vec3 col = vec3(1.0);
    float tmin = 10000.0;         // 最小の交差距離(レイの最近接交差のt値(レイがどれだけ進んだか))
    //大きい初期値を設定しているのは、レイが何も当たらない場合を示すため
    vec3  nor = vec3(0.0);        //法線
    float occ = 1.0;              //影の係数
    vec3  pos = vec3(0.0);        //交差点


    //記号整理
    //ro = r0 = レイの始点(カメラ位置)
    //rd = レイの方向
    //sc = 球の中心C(vec3(0.0,1.0,0.0))
    //t(一時変数hを利用している) = レイパラメータ(距離)
    //tmin = 現在見つかっている最小の交差距離(小さいほど手前で交差している)
    //pos = 交差点位置 = ro + t*rd
    //nor = 交差点の法線(球ならnormalize(pos - sc)、平面なら固定)

    //---------------------------
    // 地面 y=0 との交差判定
    //---------------------------
    
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
    /*
    ↑ y
   |
   |   camera (ro)
   |    *
   |     \
   |      \  rd  →  レイ方向
   |       \
   +--------+--------→ x
   y = 0  との交差  → t = (0 - ro.y) / rd.y
    */

    // ---------------------------
// 大球の交差 → fiベースで「立体小球」を物理配置
// レイは大球に当たらなくても小球に当たる可能性を評価
// ---------------------------

int   ITER       = 5;      // 反復回数（2〜6が安定）
float bumpR0     = 0.12;   // 小球半径（凸量）
float bumpRScale = 1.0;    // 反復ごとの半径スケール（1.0で固定）

vec3  accCol = vec3(0.0);  // 小球寄与の色蓄積
float accOcc = 1.0;        // 小球寄与の陰影蓄積（乗算）
float tHit    = 1e9;       // 全形状での最小ヒット距離（更新用）

// --- 大球（半径1）の交差 ---
vec3  ce = ro - sc;
float b  = dot(rd, ce);
float c  = dot(ce, ce) - 1.0;
float D  = b*b - c;

bool hitBig = false;
vec2 fi;        // 黄金螺旋インデックス＋距離
vec3 baseNor;   // 大球の法線（またはプローブ法線）
float h0 = 1e9; // 大球ヒット距離

if (D > 0.0) {
    float hEnter = -b - sqrt(D);           // 入射側
    if (hEnter > 0.0) {
        hitBig = true;
        h0     = hEnter;
        tmin   = min(tmin, hEnter);
        pos    = ro + hEnter * rd;
        baseNor= normalize(pos - sc);
        occ    = 0.5 + 0.5 * baseNor.y;
        fi     = inverseSF(baseNor);
    }
}

// --- 大球に当たらない場合でも「飛び出した小球」に当たる可能性に対応 ---
// 最接近点から単位球面の法線方向をプローブし、fiを得る
if (!hitBig) {
    float tClosest = -dot(rd, ce);               // scに対するレイの最接近パラメータ
    vec3 pClosest  = ro + rd * max(tClosest, 0.0);
    vec3 probeDir  = normalize(pClosest - sc);
    baseNor = probeDir;
    fi      = inverseSF(baseNor);
}

// --- fiに基づき、フィボナッチ点中心の「立体小球」をループで配置・交差 ---
float bumpR = bumpR0;
for (int it = 0; it < ITER; it++) {

    // 単位球面上のフィボナッチ点 q
    float id = fi.x;
    vec3 q   = fibonacciPoint(id, num);

    // 小球中心は「大球中心から q 方向へ 1.0 + bumpR」だけ外側へ
    vec3 sphC = sc + q * (1.0 + bumpR);

    // レイと小球の交差（入射側）
    vec3  ceS = ro - sphC;
    float bS  = dot(rd, ceS);
    float cS  = dot(ceS, ceS) - bumpR*bumpR;
    float DS  = bS*bS - cS;

    if (DS <= 0.0) break;                 // 交差なし → 打ち切り
    float hS = -bS - sqrt(DS);
    if (hS <= 0.0) break;                 // 後方 → 打ち切り

    // 小球ヒット
    tHit = min(tHit, hS);
    vec3 hitPos = ro + hS * rd;
    vec3 hitNor = normalize(hitPos - sphC);

    // 小球ヒット法線から「また fi を見つける」→ 次の反復用
    fi = inverseSF(hitNor);

    // 小球の色（既存のトーンに合わせて）
    vec3 iterCol = 0.5 + 0.5 * sin(hash1(fi.x * 13.0) * 3.0 + 1.0 + vec3(0.0, 1.0, 1.0));
    iterCol *= smoothstep(0.02, 0.03, fi.y);
    iterCol *= mix(1.0, 1.0 - smoothstep(0.12, 0.125, fi.y), smoothstep(-0.9, 0.1, sin(u_time)));
    iterCol *= 1.0 + 0.1 * sin(100.0 * radius * fi.y);

    // 陰影（高さで明るく）
    float iterOcc = 0.5 + 0.5 * hitNor.y;

    // 蓄積（加算＋乗算）
    accCol += iterCol * iterOcc;
    accOcc *= (0.85 + 0.15 * iterOcc);

    // 半径スケール（演出用）
    bumpR *= bumpRScale;
}

// --- 表示色統合：大球ベース＋小球凸凹寄与 ---
if (hitBig) {
    // 大球ベース色（控えめ）
    vec3 baseCol = vec3(0.9, 0.95, 1.0) * (0.5 + 0.5 * baseNor.y);
    col = baseCol + accCol;
} else {
    // 背景上に小球のみ見えるケース
    col = vec3(0.95) + accCol;
}

// 影係数の統合
occ = min(occ, accOcc);

// 最前面ヒット距離へ更新（小球が手前なら優先）
if (tHit < tmin) {
    tmin = tHit;
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
