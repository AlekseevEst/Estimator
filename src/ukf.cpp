// Ансцентный фильтр. Обобщенное ансцентное преобразование

#include "ukf.h"
#include "utils.h"

template <class M> 
unscent_filter<M>::unscent_filter()
{
    dispRgn_R = 1.0;     // дисперс. в рад.
    dispAz_R_rad = 1e-4; //
    dispUm_R_rad = 1e-4; //
    R_sph_rad.resize(3,3);
    R_sph_rad << (dispRgn_R), 0.0, 0.0, // известная ковариационная матрица ошибок измерении (СКО измерении).// ОШИБКИ ДОЛЖНЫ БЫТЬ ИЗВЕСТНЫМИ(ОШИБКИ ДАТЧИКОВ)
        0.0, (dispAz_R_rad), 0.0,
        0.0, 0.0, (dispUm_R_rad);
    
    R_sph_deg << (dispRgn_R), 0.0, 0.0, // известная ковариационная матрица ошибок измерении (СКО измерении).// ОШИБКИ ДОЛЖНЫ БЫТЬ ИЗВЕСТНЫМИ(ОШИБКИ ДАТЧИКОВ)
        0.0, (dispAz_R_rad * (180 / M_PI)), 0.0,
        0.0, 0.0, (dispUm_R_rad * (180 / M_PI));

    sa = 0.5;
}
template <class M> 
unscent_filter<M>::~unscent_filter()
{
}

template <class M> 
M unscent_filter<M>::predictUkf(const M X) // во время инициализации КФ, начальный вектор состояния = первому измерению.
{   
    // cout << " \nPREDICT FUNCTION\n ";
    // cout << "\n X = " << X << endl;

    double kappa = 1.0;
    double T = 6.0;
    Cv structCv;
    M F(6, 6); // матрица процесса (state transition matrix). (F*x(k-1)) — это как раз модель эволюции процесса.
    M G(6, 3); // матрица пересчета случайного воздействия в пространство траекторных параметров
    M Q (3,3);       // матрица ковариационных ошибок (шумовое влияние)

    if (P.isZero())
    {
        double r_meas = sqrt(pow(X(0, 0), 2) + pow(X(2, 0), 2) + pow(X(4, 0), 2));
        double az_meas =  atan2(X(2, 0), X(0, 0)) * (180/M_PI);
        double um_meas = atan2(X(4, 0), sqrt(pow(X(0, 0), 2) + pow(X(2, 0), 2))) * (180/M_PI);
        Utils<> func;
        P = func.do_cart_P(func.sph2cartcov(R_sph_deg, r_meas, az_meas, um_meas));
    }

    F << 1.0, T, 0.0, 0.0, 0.0, 0.0,  // первая строка - уравнение движения (x)
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, // вторая строка - изменение скорости (vx)
        0.0, 0.0, 1.0, T, 0.0, 0.0,   // третья строка - уравнение движения (y)
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, // четвертая строка- изменение скорости (vy)
        0.0, 0.0, 0.0, 0.0, 1.0, T,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    G << (T * T) / 2.0, 0.0, 0.0, // матрица управлении для учета внешних воздействий (xyz)
        T, 0.0, 0.0,
        0.0, (T * T) / 2.0, 0.0,
        0.0, T, 0.0,
        0.0, 0.0, (T * T) / 2.0,
        0.0, 0.0, T;

    Q << sa , 0.0, 0.0, // среднеквадратическое отклонение
        0.0, sa, 0.0,
        0.0, 0.0, sa;

    //-----------ВЗЯТИЕ МАТРИЧНОГО КОРНЯ-----------------

    Eigen::LLT<M> lltofP(P);
    if (lltofP.info() != Eigen::Success)
    {
        cout << " cholesky decomposition ERROR "; // необходима проверка матриц. чтобы не было вырождения
    }
    double n = P.cols(); // Количество сигма-векторов (столбцов в U) n равно количеству столбцов в матрице ковариации ошибок
    M U;
    M L = lltofP.matrixL();
    U = sqrt(n + kappa) * L; // Масшатбирующий коэффициент умноженный на Матричный корень
    cout << "\nU = " << U << endl;

    // Создаем словарь Xu сигма-векторов
    std::map<double, M> Xu;

    // Первый компонент
    Xu[0] = X; // в качестве первого сигма вектора берется текущий вектор состояния.
    // cout << "\nXu[0] = " << Xu[0] << endl;

    // Второй компонент. В качестве n/2 берется сумма среднего и некоторого отклонения U.col(i)
    for (double i = 0; i < n; i++)
    {
        Xu[i + 1] = X + U.col(i);
    }
    // Третий компонент. В качестве n/2 берется разность среднего и некоторого отклонения U.col(i)
    for (double i = 0; i < n; i++)
    {
        Xu[i + n + 1] = X - U.col(i);
    }

    //------------РАСЧЕТ ВЕСОВ ВЕКТОРОВ--------------
    // Первый компонент
    w[0] = kappa / (n + kappa);
    // cout << "\nw[0] = " << w[0] << endl;

    // Второй компонент
    for (double i = 0; i < 2 * n; i++)
    {
        w[i + 1] = 1.0 / (2.0 * (n + kappa));
        // cout << "\nw[ ] = " << w[i + 1] << endl;
    }
    // cout << "\nw[ ] = " << w[i+1] << endl;
    //------------ЭКСТРАПОЛЯЦИЯ СИГМА-ВЕКТОРОВ-----------------

    for (double i = 0; i < Xu.size(); i++)
    {
        Xue[i] = F * Xu[i]; // нелинейное преобразование сигма-векторов
        // cout << "\nXue[ ] = " << Xue[i] << endl;
    }
    // cout << "\nXue[0] = " << Xue[1] << endl;
    //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ----------

    structCv.Xe = M::Zero(X.rows(), X.cols());

    for (double i = 0; i < Xue.size(); i++)
    {
        structCv.Xe = structCv.Xe + w[i] * Xue[i]; // получение оценки экстр. вектора состояния
    }
    // cout << "\nstructCv.Xe = " << structCv.Xe << endl;
    //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ

    structCv.Pe = M::Zero(P.rows(), P.cols());

    for (double i = 0; i < Xue.size(); i++)
    {
        M dX;
        dX = Xue[i] - structCv.Xe;
        structCv.Pe = structCv.Pe + w[i] * dX * dX.transpose();
    }

    structCv.Pe = structCv.Pe + ((G * Q) * G.transpose());

    // cout << "\nstructCv.Pe = " << structCv.Pe << endl;

    //--------------------------------------------------------------------------------------------------
    Hcv structHcv;

    //----------ЭКСТРАПОЛИРОВАНЫЕ СИГМА-ВЕКТОРА ИЗМЕРЕНИЙ ПО НЕЛИНЕЙНЫМ ФУНКЦИЯМ------------------

    std::map<double, M> Zue;

    for (double i = 0; i < Xue.size(); i++)
    {
        M zTmp;
        zTmp.resize(3, 1);
        zTmp << sqrt(pow(Xue[i](0, 0), 2) + pow(Xue[i](2, 0), 2) + pow(Xue[i](4, 0), 2)), atan2(Xue[i](2, 0), Xue[i](0, 0)), atan2(Xue[i](4, 0), sqrt(pow(Xue[i](0, 0), 2) + pow(Xue[i](2, 0), 2)));
        Zue[i] = zTmp;
    }

    //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА ИЗМЕРЕНИЙ (ОТМЕТКИ) ------------------
    structHcv.Ze = M::Zero(3, 1);
    for (double i = 0; i < Zue.size(); i++)
    {
        structHcv.Ze = structHcv.Ze + (w[i] * Zue[i]);
    }
    // cout << "\nЭкстраполированная отметка Ze = " << structHcv.Ze << endl;

    //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА ИЗМЕРЕНИЙ
    M Pzz = M::Zero(3, 3);
    for (double i = 0; i < Zue.size(); i++)
    {
        M v;
        v = Zue[i] - structHcv.Ze;
        Pzz = Pzz + w[i] * v * v.transpose();
    }
    // cout << "\nPzz = " << Pzz << endl;
    structHcv.Se = Pzz + R_sph_rad;
    // cout << "\nstructHcv.Se = " << structHcv.Se << endl;

    //----------------------------------------------------------------------
    M Pxz = M::Zero(6, 3);

    for (double i = 0; i < Zue.size(); i++)
    {
        M dX = Xue[i] - structCv.Xe;
        M v = Zue[i] - structHcv.Ze;
        Pxz = Pxz + w[i] * dX * v.transpose();
    }

    M K = Pxz * structHcv.Se.inverse();

    predictStruct.Xe = structCv.Xe;
    predictStruct.Pe = structCv.Pe;
    predictStruct.Ze = structHcv.Ze;
    predictStruct.Se = structHcv.Se;
    predictStruct.K = K;
    M X_p(6,1);
    X_p = predictStruct.Xe;
    return X_p;
}

template <class M> 
M unscent_filter<M>::unscent_filter::correctUkf(const M Z)
{   
    // cout << " \nCORRECT FUNCTION\n ";
    // cout <<"predictStruct.Xe  = "<<predictStruct.Xe << endl;
    // cout <<"predictStruct.K = "<<predictStruct.K << endl;
    // cout <<"Z "<<Z << endl;
    // cout <<"predictStruct.Ze ="<<predictStruct.Ze << endl;

    correctStruct.X = predictStruct.Xe + predictStruct.K * (Z - predictStruct.Ze); // Z - наблюдаемые измерения(полученные),
    // cout << "\npredictStruct.Xe = \n" << predictStruct.Xe << "\n + " << " \n K*(Z-Ze)\n"<< predictStruct.K * (Z - predictStruct.Ze);
    // cout << "\n K =\n"
    //      << predictStruct.K;
    // cout << "\n Z = \n"<<Z;
    // cout << "\n Ze = \n"<<predictStruct.Ze;
    correctStruct.P = predictStruct.Pe - (predictStruct.K * predictStruct.Se) * predictStruct.K.transpose(); // по формуле
    P = correctStruct.P;
    // cout << "\nСкорректированный вектор состояния correct_struct.X = \n"
    //      << correctStruct.X;
    // cout << "\nСкорректированная матрица ошибок correct_struct.P = \n " << correct_struct.P;

    M X_c(6,1);
    X_c = correctStruct.X;
    return X_c;
}

template class unscent_filter<MatrixXd>;