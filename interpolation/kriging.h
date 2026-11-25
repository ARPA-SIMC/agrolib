/*!
 * \file kriging.h
 * \brief Interpolatore di Kriging ordinario (2D)
 *
 * Basato sull'implementazione originale di Chao-yi Lang (luglio 1995),
 * riscritto con:
 *  - variogramma con gamma(0) = 0 (il nugget agisce solo per h > 0);
 *  - risoluzione del sistema con fattorizzazione LU a pivoting parziale
 *    (la matrice di kriging ordinario ha diagonale nulla e non e' definita positiva);
 *  - varianza di kriging comprensiva del moltiplicatore di Lagrange;
 *  - controllo degli errori, nessuno stato globale nella classe, memoria gestita da std::vector;
 *  - intorno locale opzionale (k stazioni piu' vicine) per insiemi di dati grandi.
 *
 * Coordinate: le distanze sono euclidee, quindi x e y devono essere coordinate
 * cartesiane proiettate (es. UTM in metri), non latitudine/longitudine in gradi.
 *
 * Thread-safety: OrdinaryKriging::estimate() e' const e non modifica lo stato,
 * quindi piu' thread possono stimare in parallelo sulla stessa istanza.
 * Le funzioni legacy (krigingXxx) usano un'istanza globale e NON sono thread-safe.
 */

#ifndef KRIGING_H
#define KRIGING_H

#include <vector>

namespace kriging {

/*! Modelli di variogramma (i valori numerici coincidono con il vecchio parametro "mode") */
enum class Model : short
{
    Spherical   = 1,
    Exponential = 2,
    Gaussian    = 3,
    Linear      = 4
};

/*! Esito delle operazioni */
enum class Status
{
    Ok,
    NotInitialized,     //!< estimate() chiamata prima di init()
    InvalidInput,       //!< puntatori nulli, numero di stazioni < 1, valori non finiti
    InvalidParameters,  //!< parametri del variogramma non validi
    DuplicateStations,  //!< due stazioni con coordinate identiche (unirle prima)
    SingularMatrix,     //!< sistema di kriging singolare o numericamente degenere
    NegativeVariance    //!< varianza di kriging negativa oltre la tolleranza (problema numerico)
};

/*! Messaggio leggibile per uno Status */
const char* statusMessage(Status status);

/*!
 * \brief Parametri del variogramma.
 *
 * Per i modelli limitati (sferico, esponenziale, gaussiano) "sill" e' il sill TOTALE
 * (nugget compreso) e "range" e' il range pratico: a h = range il variogramma vale
 * circa il 95% del contributo strutturale (sferico: esattamente il sill).
 * Il modello lineare usa solo nugget e slope (range e sill sono ignorati).
 */
struct Variogram
{
    Model  model  = Model::Spherical;
    double range  = 1.0;
    double nugget = 0.0;
    double sill   = 1.0;
    double slope  = 1.0;

    /*! gamma(h); vale 0 per h = 0 */
    double operator()(double h) const;

    /*! controlla la coerenza dei parametri */
    Status validate() const;
};


class OrdinaryKriging
{
public:
    /*!
     * \brief Prepara l'interpolatore.
     * \param xy           coordinate delle stazioni: xy[i*2] = x, xy[i*2+1] = y
     * \param values       valori misurati, uno per stazione
     * \param count        numero di stazioni (>= 1)
     * \param variogram    modello di variogramma
     * \param maxNeighbors 0 = usa tutte le stazioni (sistema fattorizzato una sola volta);
     *                     k > 0 = per ogni punto usa solo le k stazioni piu' vicine
     *                     (consigliato quando le stazioni sono molte centinaia o piu')
     * In caso di errore l'oggetto resta non inizializzato.
     */
    Status init(const double* xy, const double* values, int count,
                const Variogram& variogram, int maxNeighbors = 0);

    /*!
     * \brief Stima il valore nel punto (x, y).
     * \param value  valore stimato (NaN in caso di errore)
     * \param stdDev se non nullo, riceve la deviazione standard di kriging (NaN in caso di errore)
     */
    Status estimate(double x, double y, double& value, double* stdDev = nullptr) const;

    /*! Libera la memoria e torna allo stato non inizializzato */
    void clear();

    bool isInitialized() const { return _ready; }
    int  size() const { return static_cast<int>(_val.size()); }
    const Variogram& variogram() const { return _variogram; }

private:
    double gammaBetween(int i, int j) const;
    void   buildMatrix(const std::vector<int>& stations, std::vector<double>& A) const;

    std::vector<double> _x, _y, _val;
    Variogram _variogram;
    int  _maxNeighbors = 0;
    bool _ready = false;

    /*! fattorizzazione LU del sistema globale (solo se _maxNeighbors == 0) */
    std::vector<double> _lu;
    std::vector<int>    _piv;
};

} // namespace kriging


/* ------------------------------------------------------------------------
 *  API legacy, compatibile con la versione precedente.
 *  Usa un'istanza globale: non thread-safe, un solo modello alla volta.
 * ------------------------------------------------------------------------ */

/*! myMode: 1 sferico, 2 esponenziale, 3 gaussiano, 4 lineare (altri valori = lineare, come prima) */
bool krigingVariogram(const double *myPos, const double *myVal, int nrItems, short myMode,
                      double myRange, double myNugget, double mySill, double mySlope);

/*! calcola stima e deviazione standard nel punto; false in caso di errore */
bool krigingSetWeight(double x_p, double y_p);

/*! valore stimato dall'ultima krigingSetWeight (NaN se fallita) */
double krigingResult();

/*! deviazione standard di kriging dall'ultima krigingSetWeight (NaN se fallita) */
double krigingStdDev();

void krigingFreeMemory();

#endif // KRIGING_H
