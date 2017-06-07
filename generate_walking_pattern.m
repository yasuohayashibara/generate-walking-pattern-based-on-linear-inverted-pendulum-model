%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;										% 変数の初期化
											% ★★★足の着地時間(s)，x方向の位置(m)，y方向の位置(m)
foot = [0 0 0; 0.8 0.0 0.2; 1.6 0.3 0.0; 2.4 0.6 0.2; 3.2 0.9 0.0; 4.0 0.9 0.2; 100 0.9 0.2];
calculate_period = 4.0;						% 歩行パターンを生成する期間(s)
dt = 0.01;									% サンプリングタイム (s)
zh = 0.8;									% 重心位置 (m)
g  = 9.8;									% 重力加速度 (m/s^2)

t = 0:dt:calculate_period;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 1;										% 目標ZMP波形の生成
n = 1;										% 足の着地時間・位置 foot → 目標ZMPパターン prefx, prefy
for tt = 0:dt:calculate_period+1
	if (tt == foot(n,1))
		prefx(i) = foot(n,2);
		prefy(i) = foot(n,3);
		n = n + 1;
	else
		prefx(i) = prefx(i-1);
		prefy(i) = prefy(i-1);
	endif
	i = i + 1;
endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = n = 1;
px = py = 0;
x = y = xi = yi = 0;
xd = yd = xdi = ydi = 0;
a = 10;
b = 1;
Tc = sqrt(zh/g);
t0 = 0;

for tt = t
	x = (xi - px)*cosh((tt - t0)/Tc)+Tc*xdi*sinh((tt - t0)/Tc)+px;
	y = (yi - py)*cosh((tt - t0)/Tc)+Tc*ydi*sinh((tt - t0)/Tc)+py;
	xd = (xi - px)/Tc*sinh((tt - t0)/Tc)+xdi*cosh((tt - t0)/Tc);
	yd = (yi - py)/Tc*sinh((tt - t0)/Tc)+ydi*cosh((tt - t0)/Tc);

	if (tt == foot(n,1))
		t0 = foot(n,1);
		xi = x;
		yi = y;
		xdi = xd;
		ydi = yd;
		px = foot(n,2);
		py = foot(n,3);
		Tsup = foot(n+1,1)-foot(n,1);
		C = cosh(Tsup/Tc);
		S = sinh(Tsup/Tc);
		xb = (foot(n+1,2)-foot(n,2))/2;
		yb = (foot(n+1,3)-foot(n,3))/2;
		vxb = (C+1)/(Tc*S)*xb;
		vyb = (C-1)/(Tc*S)*yb;
		x_des = px + xb;
		y_des = py + yb;
		xd_des = vxb;
		yd_des = vyb;
		D = a*(C-1)^2+b*(S/Tc)^2;
		px = -a*(C-1)/D*(x_des-C*xi-Tc*S*xdi)-b*S/(Tc*D)*(xd_des-S/Tc*xi-C*xdi);
		py = -a*(C-1)/D*(y_des-C*yi-Tc*S*ydi)-b*S/(Tc*D)*(yd_des-S/Tc*yi-C*ydi);
		n = n + 1;
	endif

	x0(i) = x;				% 重心位置の描画用
	y0(i) = y;
	x1(i) = prefx(i);		% 目標ZMP位置の描画用
	y1(i) = prefy(i);
	x2(i) = px;				% 実際のZMP位置の描画用
	y2(i) = py;
	i = i + 1;
endfor
subplot(2,1,1);				% グラフの描画
plot(x0, y0, "o", x1, y1, x2, y2);
subplot(2,1,2);
plot(t, x0, t, x1, t, x2, t, y0, t, y1, t, y2);
