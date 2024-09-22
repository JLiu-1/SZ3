#ifndef SZ3_QOI_INFO
#define SZ3_QOI_INFO

#include "QoI.hpp"
#include "XSquare.hpp"
#include "XCubic.hpp"
#include "XSin.hpp"
#include "XExp.hpp"
#include "XLin.hpp"
#include "LogX.hpp"
#include "XLogX.hpp"
#include "RegionalAverage.hpp"
#include "RegionalAverageOfSquare.hpp"
#include "Isoline.hpp"
#include "MultiQoIs.hpp"
#include "FX.hpp"
#include "FX_P.hpp"
#include <vector>

namespace SZ {

    template<class T, SZ::uint N>
    std::shared_ptr<concepts::QoIInterface<T, N>> GetQOI(const Config &conf){
        switch(conf.qoi){
            case 1:
                return std::make_shared<SZ::QoI_X_Square<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 2:
                return std::make_shared<SZ::QoI_Log_X<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 3:{
                if(!conf.lorenzo && !conf.lorenzo2) return std::make_shared<SZ::QoI_RegionalAverageOfSquareInterp<T, N>>(conf.qoiEB, conf.absErrorBound, conf.qoiRegionSize, conf.dims);
                else return std::make_shared<SZ::QoI_RegionalAverageOfSquare<T, N>>(conf.qoiEB, conf.absErrorBound);
                // return std::make_shared<SZ::QoI_RegionalAverage<T, N>>(conf.qoiEB, conf.absErrorBound);
            }
            case 4:{
            	std::vector<T> values;
            	for(int i=0; i<conf.isovalues.size(); i++){
            		values.push_back(conf.isovalues[i]);
            	}
                return std::make_shared<SZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound);            	
            }
            case 5:{
            	// x^2 + log x
            	std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
            	qois.push_back(std::make_shared<SZ::QoI_X_Square<T, N>>(conf.qoiEBs[0], conf.absErrorBound));
            	qois.push_back(std::make_shared<SZ::QoI_Log_X<T, N>>(conf.qoiEBs[1], conf.absErrorBound));
                return std::make_shared<SZ::QoI_MultiQoIs<T, N>>(qois);            	
            }
            case 6:{
            	// x^2 + isoline
            	std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
            	qois.push_back(std::make_shared<SZ::QoI_X_Square<T, N>>(conf.qoiEBs[0], conf.absErrorBound));
                std::vector<T> values;
                for(int i=0; i<conf.isovalues.size(); i++){
                    values.push_back(conf.isovalues[i]);
                }
                qois.push_back(std::make_shared<SZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound));
                return std::make_shared<SZ::QoI_MultiQoIs<T, N>>(qois);            	
            }
            case 7:{
                // log x + isoline
                std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
                qois.push_back(std::make_shared<SZ::QoI_Log_X<T, N>>(conf.qoiEBs[0], conf.absErrorBound));
                std::vector<T> values;
                for(int i=0; i<conf.isovalues.size(); i++){
                    values.push_back(conf.isovalues[i]);
                }
                qois.push_back(std::make_shared<SZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound));
                return std::make_shared<SZ::QoI_MultiQoIs<T, N>>(qois);             
            }
            case 8:{
            	// x^2 + log x + isoline
            	std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
            	qois.push_back(std::make_shared<SZ::QoI_X_Square<T, N>>(conf.qoiEBs[0], conf.absErrorBound));
                qois.push_back(std::make_shared<SZ::QoI_Log_X<T, N>>(conf.qoiEBs[1], conf.absErrorBound));
            	std::vector<T> values;
            	for(int i=0; i<conf.isovalues.size(); i++){
            		values.push_back(conf.isovalues[i]);
            	}
                qois.push_back(std::make_shared<SZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound));
                return std::make_shared<SZ::QoI_MultiQoIs<T, N>>(qois);            	
            }
            case 9:
                return std::make_shared<SZ::QoI_X_Cubic<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 10:
                return std::make_shared<SZ::QoI_X_Sin<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 11:
                return std::make_shared<SZ::QoI_X_Lin<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 12:
                return std::make_shared<SZ::QoI_X_Exp<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 13:
                return std::make_shared<SZ::QoI_XLog_X<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 14:
                return std::make_shared<SZ::QoI_FX<T, N>>(conf.qoiEB, conf.absErrorBound, conf.qoi_string);
            case 15:
                return std::make_shared<SZ::QoI_FX_P<T, N>>(conf.qoiEB, conf.absErrorBound, conf.qoi_string, conf.qoi_string_2, conf.threshold);
        }
        return NULL;
    }

}
#endif