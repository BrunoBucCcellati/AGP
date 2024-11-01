#include "Windows.h"
#include <queue>
#include <cliext/utility>
enum List { Top, Dawn, Left, Right };
ref class PeanoCurve_2D
{
private:
    cliext::pair<double, double>^ x1x2;
    PeanoCurve_2D^ DawnLeft;
    PeanoCurve_2D^ DawnRight;
    PeanoCurve_2D^ TopLeft;
    PeanoCurve_2D^ TopRight;
    List Type;
    double a, b, c, d;
    unsigned short razvertka;
public:
    PeanoCurve_2D(unsigned short _razvertka, List _Type, double _a, double _b, double _c, double _d)
    {
        x1x2 = gcnew cliext::pair<double, double>(0.5 * (_a + _b), 0.5 * (_c + _d)), Type = _Type, a = _a, b = _b, c = _c, d = _d, razvertka = _razvertka;
        if (_razvertka-- != 0)
        {
            if (Type == Top)
            {
                DawnLeft = gcnew PeanoCurve_2D(_razvertka, Right, _a, 0.5 * (_a + _b), _c, 0.5 * (_c + _d));
                DawnRight = gcnew PeanoCurve_2D(_razvertka, Left, 0.5 * (_a + _b), _b, _c, 0.5 * (_c + _d));
                TopLeft = gcnew PeanoCurve_2D(_razvertka, Top, _a, 0.5 * (_a + _b), 0.5 * (_c + _d), _d);
                TopRight = gcnew PeanoCurve_2D(_razvertka, Top, 0.5 * (_a + _b), _b, 0.5 * (_c + _d), _d);
            }
            if (Type == Dawn)
            {
                DawnLeft = gcnew PeanoCurve_2D(_razvertka, Dawn, _a, 0.5 * (_a + _b), _c, 0.5 * (_c + _d));
                DawnRight = gcnew PeanoCurve_2D(_razvertka, Dawn, 0.5 * (_a + _b), _b, _c, 0.5 * (_c + _d));
                TopLeft = gcnew PeanoCurve_2D(_razvertka, Right, _a, 0.5 * (_a + _b), 0.5 * (_c + _d), _d);
                TopRight = gcnew PeanoCurve_2D(_razvertka, Left, 0.5 * (_a + _b), _b, 0.5 * (_c + _d), _d);
            }
            if (Type == Left)
            {
                DawnLeft = gcnew PeanoCurve_2D(_razvertka, Left, _a, 0.5 * (_a + _b), _c, 0.5 * (_c + _d));
                DawnRight = gcnew PeanoCurve_2D(_razvertka, Top, 0.5 * (_a + _b), _b, _c, 0.5 * (_c + _d));
                TopLeft = gcnew PeanoCurve_2D(_razvertka, Left, _a, 0.5 * (_a + _b), 0.5 * (_c + _d), _d);
                TopRight = gcnew PeanoCurve_2D(_razvertka, Dawn, 0.5 * (_a + _b), _b, 0.5 * (_c + _d), _d);
            }
            if (Type == Right)
            {
                DawnLeft = gcnew PeanoCurve_2D(_razvertka, Top, _a, 0.5 * (_a + _b), _c, 0.5 * (_c + _d));
                DawnRight = gcnew PeanoCurve_2D(_razvertka, Right, 0.5 * (_a + _b), _b, _c, 0.5 * (_c + _d));
                TopLeft = gcnew PeanoCurve_2D(_razvertka, Dawn, _a, 0.5 * (_a + _b), 0.5 * (_c + _d), _d);
                TopRight = gcnew PeanoCurve_2D(_razvertka, Right, 0.5 * (_a + _b), _b, 0.5 * (_c + _d), _d);
            }
        }
    }
    cliext::pair<double, double> HitTest_2D(double x)
    {
        PeanoCurve_2D^ tmp = this;
        PeanoCurve_2D^ Curr = this;
        unsigned short num, i, _razvertka;
        do
        {
            _razvertka = tmp->razvertka;
            i = 0;
            while (i != _razvertka)
            {
                num = x * (1 << ++i + i) / (tmp->b - tmp->a);
                if (num == 0)
                {
                    if (Curr->Type == Top || Curr->Type == Right)
                    {
                        Curr = Curr->DawnLeft;
                    }
                    else
                    {
                        Curr = Curr->TopRight;
                    }
                }
                if (num == 1)
                {
                    if (Curr->Type == Top || Curr->Type == Left)
                    {
                        Curr = Curr->TopLeft;
                    }
                    else
                    {
                        Curr = Curr->DawnRight;
                    }
                }
                if (num == 2)
                {
                    if (Curr->Type == Top || Curr->Type == Right)
                    {
                        Curr = Curr->TopRight;
                    }
                    else
                    {
                        Curr = Curr->DawnLeft;
                    }
                }
                if (num == 3)
                {
                    if (Curr->Type == Top || Curr->Type == Left)
                    {
                        Curr = Curr->DawnRight;
                    }
                    else
                    {
                        Curr = Curr->TopLeft;
                    }
                }
                x -= num * (tmp->b - tmp->a) * ldexp(1, -i - i);
            }
            tmp = Curr = gcnew PeanoCurve_2D(_razvertka >> 1, Curr->Type, Curr->a, Curr->b, Curr->c, Curr->d);
        } while (_razvertka != 0);
        return cliext::pair<double, double>(Curr->x1x2);
    }
    double FindX_2D(cliext::pair<double, double> _x1x2)
    {
        PeanoCurve_2D^ tmp = this;
        PeanoCurve_2D^ Curr = this;
        unsigned short _razvertka;
        double x, x1, x2;
        do
        {
            _razvertka = tmp->razvertka;
            x = tmp->a;
            while (_razvertka != 0)
            {
                x1 = Curr->x1x2->first, x2 = Curr->x1x2->second;
                if (_x1x2.first > x1 && _x1x2.second > x2)
                {
                    if (Curr->Type == Top || Curr->Type == Right)
                    {
                        x += (tmp->b - tmp->a) * 0.5 * ldexp(1, _razvertka + _razvertka - tmp->razvertka - tmp->razvertka);
                    }
                    _razvertka--;
                    Curr = Curr->TopRight;
                }
                if (_x1x2.first < x1 && _x1x2.second > x2)
                {
                    if (Curr->Type == Top || Curr->Type == Left)
                    {
                        x += (tmp->b - tmp->a) * 0.25 * ldexp(1, _razvertka + _razvertka-- - tmp->razvertka - tmp->razvertka);
                    }
                    else
                    {
                        x += (tmp->b - tmp->a) * 0.75 * ldexp(1, _razvertka + _razvertka-- - tmp->razvertka - tmp->razvertka);
                    }
                    Curr = Curr->TopLeft;
                }
                if (_x1x2.first < x1 && _x1x2.second < x2)
                {
                    if (Curr->Type == Dawn || Curr->Type == Left)
                    {
                        x += (tmp->b - tmp->a) * 0.5 * ldexp(1, _razvertka + _razvertka - tmp->razvertka - tmp->razvertka);
                    }
                    _razvertka--;
                    Curr = Curr->DawnLeft;
                }
                if (_x1x2.first > x1 && _x1x2.second < x2)
                {
                    if (Curr->Type == Top || Curr->Type == Left)
                    {
                        x += (tmp->b - tmp->a) * 0.75 * ldexp(1, _razvertka + _razvertka-- - tmp->razvertka - tmp->razvertka);
                    }
                    else
                    {
                        x += (tmp->b - tmp->a) * 0.25 * ldexp(1, _razvertka + _razvertka-- - tmp->razvertka - tmp->razvertka);
                    }
                    Curr = Curr->DawnRight;
                }
            }
            tmp = Curr = gcnew PeanoCurve_2D(tmp->razvertka >> 1, Curr->Type, Curr->a, Curr->b, Curr->c, Curr->d);
        } while (_razvertka != 0);
        return x;
    }
};
double Sign(double Value)
{
    if (Value == 0.)
    {
        return 0;
    }
    if (Value > 0.)
    {
        return 1;
    }
    return -1;
}
double Shag(double _m, double x1, double x2, double y1, double y2, unsigned short _N)
{
    return 0.5 * (x1 + x2) - Sign(y2 - y1) * pow(0.5 * abs(y2 - y1) / _m, double(_N));
}
class Interval
{
private:
    std::pair<double, double> start, end;
    double M, R;
public:
    Interval(std::pair<double, double> _start = std::pair<double, double>(), std::pair<double, double> _end = std::pair<double, double>(), unsigned short _N = unsigned short())
    {
        start = _start, end = _end, M = abs(end.second - start.second) / pow(abs(end.first - start.first), 1 / double(_N));
    }
    void ChangeCharacteristic(double _m, unsigned short _N)
    {
        R = pow(end.first - start.first, 1 / double(_N)) + (end.second - start.second) * (end.second - start.second) / (_m * _m * pow(end.first - start.first, 1 / double(_N))) - 2 * (end.second + start.second) / _m;
    }
    double GetCharacteristic()
    {
        return R;
    }
    double GetM()
    {
        return M;
    }
    void SetCharacteristic()
    {
        R = DBL_MAX;
    }
    std::pair<double, double> GetEnd()
    {
        return end;
    }
    std::pair<double, double> GetStart()
    {
        return start;
    }
};
class Compare
{
public:
    bool operator()(Interval below, Interval above)
    {
        if (below.GetCharacteristic() < above.GetCharacteristic())
        {
            return true;
        }
        if (below.GetCharacteristic() == above.GetCharacteristic() && below.GetStart().first > above.GetStart().first)
        {
            return true;
        }
        return false;
    }
};
typedef std::priority_queue<Interval, std::vector<Interval>, Compare> Mypriority_queue;
std::deque<double> Base_LNA_1_2_Mer_AGP(bool mode = false, unsigned short N = 1, double b = 10, PeanoCurve_2D^ Curve = nullptr, PeanoCurve_2D^ Curve_Minus_PI_Na_Dva = nullptr, double r = 2, double epsilon = 5 * pow(10, -15), unsigned short global_local_iterations = 27, double a = 0, double c = 0, double d = 1)
{
    std::pair<double, double> start, end, start_Minus_PI_Na_Dva, end_Minus_PI_Na_Dva, x_Rmax, y_Rmax, x_Rmax_Minus_PI_Na_Dva, y_Rmax_Minus_PI_Na_Dva, pred_i_sled_shag, pred_i_sled_shag_Minus_PI_Na_Dva, promejutochnaya_tochka, promejutochnaya_tochka_Minus_PI_Na_Dva;
    Interval nachalny_otrezok, nachalny_otrezok_Minus_PI_Na_Dva, promejutochny_otrezok, promejutochny_otrezok_Minus_PI_Na_Dva, curr, curr1, curr_Minus_PI_Na_Dva, curr1_Minus_PI_Na_Dva, current;
    double Mmax, Mmax_Minus_PI_Na_Dva, m, m_Minus_PI_Na_Dva, dmax, dmax_Minus_PI_Na_Dva, eta_shtrih, eta_shtrih_Minus_PI_Na_Dva;
    Mypriority_queue R, R_Minus_PI_Na_Dva, R1, R1_Minus_PI_Na_Dva;
    pred_i_sled_shag = std::pair<double, double>(a, b);
    std::deque<double> Extr;
    unsigned short schetchick = 0;
    if (N == 1)
    {
        HINSTANCE load_function = LoadLibrary(L"TEST_FUNC.dll"); typedef double (*sh) (double); sh ShekelFunc = (sh)GetProcAddress(load_function, "ShekelFunc"); start = std::pair<double, double>(a, ShekelFunc(a)), end = std::pair<double, double>(b, ShekelFunc(b)), nachalny_otrezok = Interval(start, end, N), Mmax = nachalny_otrezok.GetM(), m = r * Mmax, x_Rmax = std::pair<double, double>(start.first, end.first), y_Rmax = std::pair<double, double>(start.second, end.second), R.push(nachalny_otrezok);
        while (abs(pred_i_sled_shag.second - pred_i_sled_shag.first) > epsilon)
        {
            pred_i_sled_shag.first = pred_i_sled_shag.second, promejutochnaya_tochka.first = pred_i_sled_shag.second = Shag(m, x_Rmax.first, x_Rmax.second, y_Rmax.first, y_Rmax.second, N), Extr.push_back(ShekelFunc(pred_i_sled_shag.second));
            if (schetchick == 1000)
            {
                return Extr;
            }
            promejutochnaya_tochka.second = Extr.back(), promejutochny_otrezok = R.top(), curr = Interval(promejutochny_otrezok.GetStart(), promejutochnaya_tochka, N), curr1 = Interval(promejutochnaya_tochka, promejutochny_otrezok.GetEnd(), N), R.pop();
            if (mode == true && schetchick > global_local_iterations && schetchick % 2 == 0)
            {
                while (R.empty() == false)
                {
                    promejutochny_otrezok = R.top(), eta_shtrih = (std::max)((std::max)(curr.GetM(), curr1.GetM()), promejutochny_otrezok.GetM()), promejutochny_otrezok.ChangeCharacteristic(r * Mmax * (promejutochny_otrezok.GetEnd().first - promejutochny_otrezok.GetStart().first) * dmax + eta_shtrih - Mmax * (promejutochny_otrezok.GetEnd().first - promejutochny_otrezok.GetStart().first) * dmax, N), R1.push(promejutochny_otrezok), R.pop();
                }
                R = R1, R1 = Mypriority_queue();
            }
            else
            {
                if ((std::max)(curr.GetM(), curr1.GetM()) < Mmax)
                {
                    curr.ChangeCharacteristic(m, N), curr1.ChangeCharacteristic(m, N);
                }
                else
                {
                    Mmax = (std::max)(curr.GetM(), curr1.GetM()), m = r * Mmax, curr.ChangeCharacteristic(m, N), curr1.ChangeCharacteristic(m, N);
                    if (mode == true)
                    {
                        dmax = (std::max)(pow((curr.GetEnd()).first - (curr.GetStart()).first, (1 / double(N))), pow((curr1.GetEnd()).first - (curr1.GetStart()).first, (1 / double(N))));
                    }
                    while (R.empty() == false)
                    {
                        promejutochny_otrezok = R.top();
                        if (mode == true && pow((promejutochny_otrezok.GetEnd()).first - (promejutochny_otrezok.GetStart()).first, (1 / double(N))) > dmax)
                        {
                            dmax = pow((promejutochny_otrezok.GetEnd()).first - (promejutochny_otrezok.GetStart()).first, (1 / double(N)));
                        }
                        promejutochny_otrezok.ChangeCharacteristic(m, N), R1.push(promejutochny_otrezok), R.pop();
                    }
                    R = R1, R1 = Mypriority_queue();
                }
            }
            R.push(curr), R.push(curr1), promejutochny_otrezok = R.top(), x_Rmax.first = promejutochny_otrezok.GetStart().first, x_Rmax.second = promejutochny_otrezok.GetEnd().first, y_Rmax.first = promejutochny_otrezok.GetStart().second, y_Rmax.second = promejutochny_otrezok.GetEnd().second, schetchick++;
        }
        FreeLibrary(load_function);
    }
    else
    {
        HINSTANCE load_function = LoadLibrary(L"TEST_FUNC.dll"); typedef double (*grsh) (double, double); grsh GrishaginFunc = (grsh)GetProcAddress(load_function, "GrishaginFunc"); start = std::pair<double, double>(a, GrishaginFunc(a, c)), end = std::pair<double, double>(b, GrishaginFunc(b, c)), start_Minus_PI_Na_Dva = std::pair<double, double>(a, GrishaginFunc(a, c)), end_Minus_PI_Na_Dva = std::pair<double, double>(b, GrishaginFunc(a, d));
        nachalny_otrezok = Interval(start, end, N), nachalny_otrezok_Minus_PI_Na_Dva = Interval(start_Minus_PI_Na_Dva, end_Minus_PI_Na_Dva, N), Mmax = nachalny_otrezok.GetM(), m = r * Mmax, x_Rmax = std::pair<double, double>(start.first, end.first), y_Rmax = std::pair<double, double>(start.second, end.second), Mmax_Minus_PI_Na_Dva = nachalny_otrezok_Minus_PI_Na_Dva.GetM(), m_Minus_PI_Na_Dva = r * Mmax_Minus_PI_Na_Dva, R.push(nachalny_otrezok), R_Minus_PI_Na_Dva.push(nachalny_otrezok_Minus_PI_Na_Dva), x_Rmax_Minus_PI_Na_Dva = std::pair<double, double>(start_Minus_PI_Na_Dva.first, end_Minus_PI_Na_Dva.first), y_Rmax_Minus_PI_Na_Dva = std::pair<double, double>(start_Minus_PI_Na_Dva.second, end_Minus_PI_Na_Dva.second), pred_i_sled_shag_Minus_PI_Na_Dva = std::pair<double, double>(a, b);
        while (abs(pred_i_sled_shag.second - pred_i_sled_shag.first) > epsilon)
        {
            pred_i_sled_shag.first = pred_i_sled_shag.second, pred_i_sled_shag.second = Shag(m, x_Rmax.first, x_Rmax.second, y_Rmax.first, y_Rmax.second, N), pred_i_sled_shag_Minus_PI_Na_Dva.first = pred_i_sled_shag_Minus_PI_Na_Dva.second, pred_i_sled_shag_Minus_PI_Na_Dva.second = Shag(m_Minus_PI_Na_Dva, x_Rmax_Minus_PI_Na_Dva.first, x_Rmax_Minus_PI_Na_Dva.second, y_Rmax_Minus_PI_Na_Dva.first, y_Rmax_Minus_PI_Na_Dva.second, N);
            if (schetchick % 20 == 0 && schetchick != 0)
            {
                Extr.push_back((std::min)(GrishaginFunc(Curve->HitTest_2D(pred_i_sled_shag.second).first, Curve->HitTest_2D(pred_i_sled_shag.second).second), GrishaginFunc(Curve_Minus_PI_Na_Dva->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second).first, Curve_Minus_PI_Na_Dva->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second).second)));
                if (schetchick == 1000)
                {
                    return Extr;
                }
                if (Extr.back() == GrishaginFunc(Curve->HitTest_2D(pred_i_sled_shag.second).first, Curve->HitTest_2D(pred_i_sled_shag.second).second))
                {
                    pred_i_sled_shag_Minus_PI_Na_Dva.second = Curve_Minus_PI_Na_Dva->FindX_2D(Curve->HitTest_2D(pred_i_sled_shag.second));
                    while (!(promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first - epsilon < pred_i_sled_shag_Minus_PI_Na_Dva.second && pred_i_sled_shag_Minus_PI_Na_Dva.second < promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first + epsilon))
                    {
                        promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(), R_Minus_PI_Na_Dva.pop();
                        if (promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first - epsilon < pred_i_sled_shag_Minus_PI_Na_Dva.second && pred_i_sled_shag_Minus_PI_Na_Dva.second < promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first + epsilon)
                        {
                            promejutochny_otrezok_Minus_PI_Na_Dva.SetCharacteristic();
                        }
                        R1.push(promejutochny_otrezok_Minus_PI_Na_Dva);
                    }
                    while (R1.empty() == false)
                    {
                        current = R1.top(), R1.pop(), R_Minus_PI_Na_Dva.push(current);
                    }
                }
                else
                {
                    pred_i_sled_shag.second = Curve->FindX_2D(Curve_Minus_PI_Na_Dva->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second));
                    while (!(promejutochny_otrezok.GetStart().first - epsilon < pred_i_sled_shag.second && pred_i_sled_shag.second < promejutochny_otrezok.GetEnd().first + epsilon))
                    {
                        promejutochny_otrezok = R.top(), R.pop();
                        if (promejutochny_otrezok.GetStart().first - epsilon < pred_i_sled_shag.second && pred_i_sled_shag.second < promejutochny_otrezok.GetEnd().first + epsilon)
                        {
                            promejutochny_otrezok.SetCharacteristic();
                        }
                        R1.push(promejutochny_otrezok);
                    }
                    while (R1.empty() == false)
                    {
                        current = R1.top(), R1.pop(), R.push(current);
                    }
                }
                promejutochnaya_tochka.second = promejutochnaya_tochka_Minus_PI_Na_Dva.second = Extr.back();
            }
            else
            {
                promejutochnaya_tochka.second = GrishaginFunc(Curve->HitTest_2D(pred_i_sled_shag.second).first, Curve->HitTest_2D(pred_i_sled_shag.second).second); promejutochnaya_tochka_Minus_PI_Na_Dva.second = GrishaginFunc(Curve_Minus_PI_Na_Dva->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second).first, Curve_Minus_PI_Na_Dva->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second).second), promejutochny_otrezok = R.top(), promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(), Extr.push_back((std::min)(promejutochnaya_tochka.second, promejutochnaya_tochka_Minus_PI_Na_Dva.second));
            }
            promejutochnaya_tochka.first = pred_i_sled_shag.second; promejutochnaya_tochka_Minus_PI_Na_Dva.first = pred_i_sled_shag_Minus_PI_Na_Dva.second, curr = Interval(promejutochny_otrezok.GetStart(), promejutochnaya_tochka, N), curr1 = Interval(promejutochnaya_tochka, promejutochny_otrezok.GetEnd(), N), R.pop(), curr_Minus_PI_Na_Dva = Interval(promejutochny_otrezok_Minus_PI_Na_Dva.GetStart(), promejutochnaya_tochka_Minus_PI_Na_Dva, N), curr1_Minus_PI_Na_Dva = Interval(promejutochnaya_tochka_Minus_PI_Na_Dva, promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd(), N), R_Minus_PI_Na_Dva.pop();
            if (mode == true && schetchick > global_local_iterations && schetchick % 2 == 0)
            {
                while (R.empty() == false)
                {
                    promejutochny_otrezok = R.top(), eta_shtrih = (std::max)((std::max)(curr.GetM(), curr1.GetM()), promejutochny_otrezok.GetM()), promejutochny_otrezok.ChangeCharacteristic(r * Mmax * (promejutochny_otrezok.GetEnd().first - promejutochny_otrezok.GetStart().first) * dmax + eta_shtrih - Mmax * (promejutochny_otrezok.GetEnd().first - promejutochny_otrezok.GetStart().first) * dmax, N), R1.push(promejutochny_otrezok), R.pop();
                }
                R = R1, R1 = Mypriority_queue();
                while (R_Minus_PI_Na_Dva.empty() == false)
                {
                    promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(), eta_shtrih_Minus_PI_Na_Dva = (std::max)((std::max)(curr_Minus_PI_Na_Dva.GetM(), curr1_Minus_PI_Na_Dva.GetM()), promejutochny_otrezok_Minus_PI_Na_Dva.GetM()), promejutochny_otrezok_Minus_PI_Na_Dva.ChangeCharacteristic(r * Mmax_Minus_PI_Na_Dva * (promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first - promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first) * dmax_Minus_PI_Na_Dva + eta_shtrih_Minus_PI_Na_Dva - Mmax_Minus_PI_Na_Dva * (promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first - promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first) * dmax_Minus_PI_Na_Dva, N), R1_Minus_PI_Na_Dva.push(promejutochny_otrezok_Minus_PI_Na_Dva), R_Minus_PI_Na_Dva.pop();
                }
                R_Minus_PI_Na_Dva = R1_Minus_PI_Na_Dva, R1_Minus_PI_Na_Dva = Mypriority_queue();
            }
            else
            {
                if ((std::max)(curr.GetM(), curr1.GetM()) < Mmax)
                {
                    curr.ChangeCharacteristic(m, N), curr1.ChangeCharacteristic(m, N);
                }
                else
                {
                    Mmax = (std::max)(curr.GetM(), curr1.GetM()), m = r * Mmax, curr.ChangeCharacteristic(m, N), curr1.ChangeCharacteristic(m, N);
                    if (mode == true)
                    {
                        dmax = (std::max)(pow((curr.GetEnd()).first - (curr.GetStart()).first, (1 / double(N))), pow((curr1.GetEnd()).first - (curr1.GetStart()).first, (1 / double(N))));
                    }
                    while (R.empty() == false)
                    {
                        promejutochny_otrezok = R.top();
                        if (mode == true && pow((promejutochny_otrezok.GetEnd()).first - (promejutochny_otrezok.GetStart()).first, (1 / double(N))) > dmax)
                        {
                            dmax = pow((promejutochny_otrezok.GetEnd()).first - (promejutochny_otrezok.GetStart()).first, (1 / double(N)));
                        }
                        promejutochny_otrezok.ChangeCharacteristic(m, N), R1.push(promejutochny_otrezok), R.pop();
                    }
                    R = R1, R1 = Mypriority_queue();
                }
                if ((std::max)(curr_Minus_PI_Na_Dva.GetM(), curr1_Minus_PI_Na_Dva.GetM()) < Mmax_Minus_PI_Na_Dva)
                {
                    curr_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N), curr1_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N);
                }
                else
                {
                    Mmax_Minus_PI_Na_Dva = (std::max)(curr_Minus_PI_Na_Dva.GetM(), curr1_Minus_PI_Na_Dva.GetM()), m_Minus_PI_Na_Dva = r * Mmax_Minus_PI_Na_Dva, curr_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N), curr1_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N);
                    if (mode == true)
                    {
                        dmax_Minus_PI_Na_Dva = (std::max)(pow((curr_Minus_PI_Na_Dva.GetEnd()).first - (curr_Minus_PI_Na_Dva.GetStart()).first, (1 / double(N))), pow((curr1_Minus_PI_Na_Dva.GetEnd()).first - (curr1_Minus_PI_Na_Dva.GetStart()).first, (1 / double(N))));
                    }
                    while (R_Minus_PI_Na_Dva.empty() == false)
                    {
                        promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top();
                        if (mode == true && pow((promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd()).first - (promejutochny_otrezok_Minus_PI_Na_Dva.GetStart()).first, (1 / double(N))) > dmax_Minus_PI_Na_Dva)
                        {
                            dmax_Minus_PI_Na_Dva = pow((promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd()).first - (promejutochny_otrezok_Minus_PI_Na_Dva.GetStart()).first, (1 / double(N)));
                        }
                        promejutochny_otrezok_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N), R1_Minus_PI_Na_Dva.push(promejutochny_otrezok_Minus_PI_Na_Dva), R_Minus_PI_Na_Dva.pop();
                    }
                    R_Minus_PI_Na_Dva = R1_Minus_PI_Na_Dva, R1_Minus_PI_Na_Dva = Mypriority_queue();
                }
            }
            R.push(curr), R.push(curr1), promejutochny_otrezok = R.top(), x_Rmax.first = promejutochny_otrezok.GetStart().first, x_Rmax.second = promejutochny_otrezok.GetEnd().first, y_Rmax.first = promejutochny_otrezok.GetStart().second, y_Rmax.second = promejutochny_otrezok.GetEnd().second, R_Minus_PI_Na_Dva.push(curr_Minus_PI_Na_Dva), R_Minus_PI_Na_Dva.push(curr1_Minus_PI_Na_Dva), promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(), x_Rmax_Minus_PI_Na_Dva.first = promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first, x_Rmax_Minus_PI_Na_Dva.second = promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first, y_Rmax_Minus_PI_Na_Dva.first = promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().second, y_Rmax_Minus_PI_Na_Dva.second = promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().second, schetchick++;
        }
        FreeLibrary(load_function);
    }
    return Extr;
}