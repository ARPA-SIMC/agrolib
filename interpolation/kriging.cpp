/*!
 * \file kriging.cpp
 * \brief Interpolatore di Kriging ordinario (2D) - vedi kriging.h
 */

#include "kriging.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace kriging {

namespace {

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

/*!
 * \brief Fattorizzazione LU in place con pivoting parziale (PA = LU).
 * \param A   matrice n*n row-major; in uscita contiene L (sotto la diagonale, diagonale unitaria) e U
 * \param piv in uscita: piv[k] = riga scambiata con la riga k al passo k
 * \return false se la matrice e' singolare o numericamente degenere
 */
bool luDecompose(std::vector<double>& A, std::vector<int>& piv, int n)
{
    piv.assign(n, 0);

    double maxAbs = 0.;
    for (double a : A)
        maxAbs = std::max(maxAbs, std::fabs(a));
    if (! (maxAbs > 0.) || ! std::isfinite(maxAbs))
        return false;

    const double tolerance = n * std::numeric_limits<double>::epsilon() * maxAbs;

    for (int k = 0; k < n; k++)
    {
        /* ricerca del pivot sulla colonna k */
        int p = k;
        double best = std::fabs(A[k*n + k]);
        for (int i = k + 1; i < n; i++)
        {
            const double v = std::fabs(A[i*n + k]);
            if (v > best)
            {
                best = v;
                p = i;
            }
        }
        if (best <= tolerance)
            return false;

        piv[k] = p;
        if (p != k)
            for (int j = 0; j < n; j++)
                std::swap(A[k*n + j], A[p*n + j]);

        /* eliminazione */
        const double pivot = A[k*n + k];
        for (int i = k + 1; i < n; i++)
        {
            const double factor = (A[i*n + k] /= pivot);
            if (factor == 0.)
                continue;
            for (int j = k + 1; j < n; j++)
                A[i*n + j] -= factor * A[k*n + j];
        }
    }
    return true;
}

/*! \brief Risolve LU x = P b; b viene sovrascritto con la soluzione */
void luSolve(const std::vector<double>& LU, const std::vector<int>& piv, std::vector<double>& b, int n)
{
    for (int k = 0; k < n; k++)
        if (piv[k] != k)
            std::swap(b[k], b[piv[k]]);

    /* sostituzione in avanti (L ha diagonale unitaria) */
    for (int i = 1; i < n; i++)
    {
        double s = b[i];
        for (int j = 0; j < i; j++)
            s -= LU[i*n + j] * b[j];
        b[i] = s;
    }

    /* sostituzione all'indietro */
    for (int i = n - 1; i >= 0; i--)
    {
        double s = b[i];
        for (int j = i + 1; j < n; j++)
            s -= LU[i*n + j] * b[j];
        b[i] = s / LU[i*n + i];
    }
}

} // namespace


const char* statusMessage(Status status)
{
    switch (status)
    {
        case Status::Ok:                return "ok";
        case Status::NotInitialized:    return "interpolatore non inizializzato";
        case Status::InvalidInput:      return "dati di input non validi";
        case Status::InvalidParameters: return "parametri del variogramma non validi";
        case Status::DuplicateStations: return "stazioni con coordinate duplicate";
        case Status::SingularMatrix:    return "sistema di kriging singolare";
        case Status::NegativeVariance:  return "varianza di kriging negativa (instabilita' numerica)";
    }
    return "errore sconosciuto";
}


/* ---------------------------------------------------------------- Variogram */

double Variogram::operator()(double h) const
{
    /* gamma(0) = 0 per definizione: il nugget e' una discontinuita' all'origine */
    if (h <= 0.)
        return 0.;

    const double structural = sill - nugget;
    const double r = h / range;

    switch (model)
    {
        case Model::Spherical:
            return (h < range) ? nugget + structural * (1.5 * r - 0.5 * r * r * r)
                               : sill;

        /* range pratico: fattore 3 (convenzione GSLIB), gamma(range) ~ 95% */
        case Model::Exponential:
            return nugget + structural * (1. - std::exp(-3. * r));

        case Model::Gaussian:
            return nugget + structural * (1. - std::exp(-3. * r * r));

        case Model::Linear:
            return nugget + slope * h;
    }
    return NaN;
}


Status Variogram::validate() const
{
    if (! std::isfinite(nugget) || nugget < 0.)
        return Status::InvalidParameters;

    switch (model)
    {
        case Model::Linear:
            if (! std::isfinite(slope) || slope < 0.)
                return Status::InvalidParameters;
            return Status::Ok;

        case Model::Spherical:
        case Model::Exponential:
        case Model::Gaussian:
            if (! std::isfinite(range) || range <= 0.)
                return Status::InvalidParameters;
            if (! std::isfinite(sill) || sill < nugget)
                return Status::InvalidParameters;
            return Status::Ok;
    }
    return Status::InvalidParameters;
}


/* ---------------------------------------------------------- OrdinaryKriging */

double OrdinaryKriging::gammaBetween(int i, int j) const
{
    return _variogram(std::hypot(_x[i] - _x[j], _y[i] - _y[j]));
}


/*!
 * \brief Costruisce la matrice di kriging ordinario (m+1)x(m+1) sulle stazioni indicate:
 *        | gamma_ij  1 |
 *        |   1       0 |
 */
void OrdinaryKriging::buildMatrix(const std::vector<int>& stations, std::vector<double>& A) const
{
    const int m = static_cast<int>(stations.size());
    const int n = m + 1;
    A.assign(static_cast<size_t>(n) * n, 0.);

    for (int i = 0; i < m; i++)
    {
        for (int j = i + 1; j < m; j++)
            A[i*n + j] = A[j*n + i] = gammaBetween(stations[i], stations[j]);
        A[i*n + m] = 1.;
        A[m*n + i] = 1.;
    }
    /* diagonale e angolo in basso a destra restano a 0 */
}


Status OrdinaryKriging::init(const double* xy, const double* values, int count,
                             const Variogram& variogram, int maxNeighbors)
{
    clear();

    if (xy == nullptr || values == nullptr || count < 1 || maxNeighbors < 0)
        return Status::InvalidInput;

    const Status vStatus = variogram.validate();
    if (vStatus != Status::Ok)
        return vStatus;

    _x.resize(count);
    _y.resize(count);
    _val.resize(count);
    for (int i = 0; i < count; i++)
    {
        _x[i]   = xy[i*2];
        _y[i]   = xy[i*2 + 1];
        _val[i] = values[i];
        if (! std::isfinite(_x[i]) || ! std::isfinite(_y[i]) || ! std::isfinite(_val[i]))
        {
            clear();
            return Status::InvalidInput;
        }
    }

    /* stazioni duplicate: righe identiche nella matrice -> sistema singolare */
    std::vector<int> order(count);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [this](int a, int b)
              { return _x[a] < _x[b] || (_x[a] == _x[b] && _y[a] < _y[b]); });
    for (int k = 1; k < count; k++)
    {
        if (_x[order[k]] == _x[order[k-1]] && _y[order[k]] == _y[order[k-1]])
        {
            clear();
            return Status::DuplicateStations;
        }
    }

    _variogram = variogram;
    _maxNeighbors = (maxNeighbors >= count) ? 0 : maxNeighbors;

    /* sistema globale: fattorizzato una sola volta, O(n^3) */
    if (_maxNeighbors == 0)
    {
        std::vector<int> all(count);
        std::iota(all.begin(), all.end(), 0);
        buildMatrix(all, _lu);
        if (! luDecompose(_lu, _piv, count + 1))
        {
            clear();
            return Status::SingularMatrix;
        }
    }

    _ready = true;
    return Status::Ok;
}


Status OrdinaryKriging::estimate(double x, double y, double& value, double* stdDev) const
{
    value = NaN;
    if (stdDev != nullptr)
        *stdDev = NaN;

    if (! _ready)
        return Status::NotInitialized;
    if (! std::isfinite(x) || ! std::isfinite(y))
        return Status::InvalidInput;

    const int count = size();

    /* stazioni usate per questo punto */
    std::vector<int> stations(count);
    std::iota(stations.begin(), stations.end(), 0);

    if (_maxNeighbors > 0)
    {
        std::vector<double> dist2(count);
        for (int i = 0; i < count; i++)
        {
            const double dx = _x[i] - x;
            const double dy = _y[i] - y;
            dist2[i] = dx * dx + dy * dy;
        }
        std::nth_element(stations.begin(), stations.begin() + _maxNeighbors, stations.end(),
                         [&dist2](int a, int b) { return dist2[a] < dist2[b]; });
        stations.resize(_maxNeighbors);
    }

    const int m = static_cast<int>(stations.size());
    const int n = m + 1;

    /* termine noto: gamma(stazione, P) e vincolo sum(w) = 1 */
    std::vector<double> gamma0(m);
    std::vector<double> b(n);
    for (int i = 0; i < m; i++)
    {
        const int s = stations[i];
        gamma0[i] = _variogram(std::hypot(_x[s] - x, _y[s] - y));
        b[i] = gamma0[i];
    }
    b[m] = 1.;

    /* soluzione: b = [w_0 .. w_m-1, mu] */
    if (_maxNeighbors > 0)
    {
        std::vector<double> lu;
        std::vector<int> piv;
        buildMatrix(stations, lu);
        if (! luDecompose(lu, piv, n))
            return Status::SingularMatrix;
        luSolve(lu, piv, b, n);
    }
    else
    {
        luSolve(_lu, _piv, b, n);
    }

    double estimateValue = 0.;
    double variance = b[m];                 /* moltiplicatore di Lagrange mu */
    double magnitude = std::fabs(b[m]);     /* scala del problema, per la tolleranza sulla varianza */
    for (int i = 0; i < m; i++)
    {
        estimateValue += b[i] * _val[stations[i]];
        variance  += b[i] * gamma0[i];
        /* gamma0 e non w*gamma0: su una stazione w ~ (0..1..0) e i prodotti sono tutti ~0,
           la tolleranza crollerebbe al livello del rumore di arrotondamento */
        magnitude += std::fabs(gamma0[i]);
    }

    if (! std::isfinite(estimateValue) || ! std::isfinite(variance))
        return Status::SingularMatrix;

    value = estimateValue;

    /* piccoli valori negativi sono rumore di arrotondamento, valori grandi no */
    if (variance < -1e-9 * magnitude)
        return Status::NegativeVariance;

    if (stdDev != nullptr)
        *stdDev = std::sqrt(std::max(0., variance));

    return Status::Ok;
}


void OrdinaryKriging::clear()
{
    std::vector<double>().swap(_x);
    std::vector<double>().swap(_y);
    std::vector<double>().swap(_val);
    std::vector<double>().swap(_lu);
    std::vector<int>().swap(_piv);
    _maxNeighbors = 0;
    _ready = false;
}

} // namespace kriging


/* ------------------------------------------------------------------ API legacy */

namespace {
kriging::OrdinaryKriging legacyKriging;
double legacyValue  = std::numeric_limits<double>::quiet_NaN();
double legacyStdDev = std::numeric_limits<double>::quiet_NaN();
}

bool krigingVariogram(const double *myPos, const double *myVal, int nrItems, short myMode,
                      double myRange, double myNugget, double mySill, double mySlope)
{
    legacyValue = legacyStdDev = std::numeric_limits<double>::quiet_NaN();

    kriging::Variogram v;
    /* come nella versione originale, un mode sconosciuto equivale al lineare */
    v.model  = (myMode >= 1 && myMode <= 4) ? static_cast<kriging::Model>(myMode)
                                             : kriging::Model::Linear;
    v.range  = myRange;
    v.nugget = myNugget;
    v.sill   = mySill;
    v.slope  = mySlope;

    return legacyKriging.init(myPos, myVal, nrItems, v) == kriging::Status::Ok;
}

bool krigingSetWeight(double x_p, double y_p)
{
    const kriging::Status s = legacyKriging.estimate(x_p, y_p, legacyValue, &legacyStdDev);
    if (s != kriging::Status::Ok)
    {
        legacyValue = legacyStdDev = std::numeric_limits<double>::quiet_NaN();
        return false;
    }
    return true;
}

double krigingResult()
{
    return legacyValue;
}

double krigingStdDev()
{
    return legacyStdDev;
}

void krigingFreeMemory()
{
    legacyKriging.clear();
    legacyValue = legacyStdDev = std::numeric_limits<double>::quiet_NaN();
}
