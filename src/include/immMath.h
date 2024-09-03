// #pragma once
// #include "utils.h"
// #include <vector>
// template <class M>
// struct ImmMath
// {
//     M computeMixingProbability(const M &p_ij, const M &mu_i)
//     {
//         cj.resize(1, mu_i.cols());

//         for (long int j = 0; j < mu_i.cols(); j++)
//         {
//             double c = 0.;
//             for (long int i = 0; i < p_ij.cols(); i++)
//             {
//                 c +=(p_ij(i, j) * mu_i(0, i));
//             }
//             cj(0, j) = c;
//         }

//         M mix(p_ij.rows(), p_ij.cols());

//         for (long int j = 0; j < p_ij.cols(); j++)
//         {
//             for (long int i = 0; i < p_ij.rows(); i++)
//             {
//                 mix(i, j) = p_ij(i, j) * mu_i(0, i) / cj(0, j);
//             }
//         }
//         return mix;
//     }

//     std::pair<std::vector<M>, std::vector<M>> InitMixingStateAndCovariance(const M &mu_ij)
//     {
//         std::vector<M> statesOfFilters;
//         std::vector<M> covarianceOfFilters;

//         for (size_t j = 0; j < conteiner.filters.size(); ++j)
//         {
//             M mixedState = M::Zero(conteiner.filters[j]->correctStruct.X.rows(), conteiner.filters[j]->correctStruct.X.cols());
//             M mixedCovariance = M::Zero(conteiner.filters[j]->correctStruct.P.rows(), conteiner.filters[j]->correctStruct.P.cols());

//             for (size_t i = 0; i < conteiner.filters.size(); ++i)
//             {
//                 M convertedState = converter.m[{conteiner.filters[i]->getModelType(),conteiner.filters[j]->getModelType()}](conteiner.filters[i]->correctStruct.X);
//                 mixedState += mu_ij(i, j) * convertedState;
//             }
//                 statesOfFilters.push_back(mixedState);

//             for (size_t i = 0; i < conteiner.filters.size(); ++i)
//             {
//                 M convertedState = converter.m[{conteiner.filters[i]->getModelType(),conteiner.filters[j]->getModelType()}](conteiner.filters[i]->correctStruct.X);
//                 M dX = convertedState - mixedState;
//                 M convertedCovariance = converter.m[{conteiner.filters[i]->getModelType(),conteiner.filters[j]->getModelType()}](conteiner.filters[i]->correctStruct.P);
//                 mixedCovariance += mu_ij(i, j) * (convertedCovariance + dX * dX.transpose());
//             }
                        
//             covarianceOfFilters.push_back(mixedCovariance);
//         }
//         return std::make_pair(statesOfFilters, covarianceOfFilters);
//     }
// };
    
