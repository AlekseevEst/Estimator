// #include <catch2/catch.hpp>
// #include <iostream>
// #include <Eigen/Dense>
// #include <chrono>
// #include <cinttypes>
// #include <iostream>
// #include <signal.h>
// #include <thread>
// #include "channel.h"
// #include "track.h"
// #include "JVC.h"
// #include "imm.h"
// #include "initFilters.h"

// using namespace Catch::Benchmark;
// using M = Eigen::MatrixXd;
// // struct Detection {};
// struct SectorLabel {
//     double Time;
// };
// struct CloseTask {};

// /**
// * Базовый класс, обёртка для диспечера, приём и отправка кодограмм
// */
// struct CdgParser {
//     void Run() {
//         try {
//             // запускаешь диспечер (Рахманова)
//             for (;;) {
//                 // тут разбераешь кодограммы при приёме
//                 // и отправляешь их в другой поток в tracker
//             }
//         } catch(CloseTask const&) {
//             std::cout << "[CdgParser] Done " << std::endl;
//         }
//     }

//     void Done() {
//         GetSender().Send(CloseTask());
//     }

//     Channel::Sender GetSender() {
//         return incoming;
//     }

//     Channel::Sender ToTrackerGNN;

// private:

//     Channel::Receiver incoming;

// };


// template <class TypeDetection,
//           class TypeTrack>
// struct BinderRude {
//     bool operator()(const TypeTrack& track,
//                     const TypeDetection& detection) {
//         return true;
//     }
// };

// template <class TypeDetection,
//           class TypeTrack>
// struct Binder {
//     double operator()(const TypeTrack& track,
//                       const TypeDetection& detection) {
//                         double dist;
//         return dist; // считается дистанция
//     }
// };


// template <class TypeAlgorithm,
//           class TypeDetection,
//           class TypeTrack>
// struct Association {

//     void operator()(TypeTrack& track,
//                     TypeDetection& detection, //const
//                     double distance) {
//         dictDistance[{&track, &detection}] = distance;
//         dictTracks[&detection].push_back(&track); // То есть к каждой отметке привязывается трасса, а не наоборот?
//     }

//     void Step(SectorLabel const& label, std::vector<TypeTrack*>& tracksChange) {
//         //формирую класстеры, через общие измерения, алгоритм поиска в ширину, обрабатываю dictTracks
//         //проверяю про секторная метка прошла по класстерам
//         //делаю шаг
//         //удаляю класстер
//     }

//     bool checkTrackAssociation(TypeTrack &track)
//     {
//         // for (const auto &pair : dictTracks)
//         // {
//         //     const std::vector<TypeTrack *> &tracks = pair.second;

//         //     for (auto &tr : tracks)
//         //     {
//         //         if (tr == track)
//         //             return true;
//         //     }
//         // }
//         return false;
//     }

// private:
//     TypeAlgorithm algorithm; // JVC, аукционна, венгерский и т.д.

//     std::unordered_map<std::pair<TypeTrack*, TypeDetection*>, double> dictDistance;
//     std::unordered_map<TypeDetection*, std::vector<TypeTrack*>> dictTracks;

// };
    
//     template <class TypeTrack>
//     struct Cover{
//         bool operator()(TypeTrack&,  const SectorLabel&)
//         {
//             return true;
//         }
//     };
//     template <class TypeTrack>
//     struct FormerToUser{
//         void ToUser(std::vector<TypeTrack*> tracks){

//         }
//         void ToUserDeleted(std::vector<TypeTrack*> tracks){

//         };
//     };

//     template<class TypeTrack>
//     struct StorageValidTrack : public std::vector<TypeTrack> {

//         TypeTrack GetTrack()
//         {

//         }

//         void InsertTrack(TypeTrack& track)
//         {
            
//         }

//         void Delete(std::vector<TypeTrack*>)
//         {

//         }

//     };

//     template<class TypeTrack> 
//     struct StorageTrack : public std::vector<TypeTrack>{
//         TypeTrack GetTrack(){
//             return TypeTrack();

//         }
//         void Insert(std::vector<TypeTrack*>){
            
//         }  

//     };
// /*
//     Базовый класс, алгоритм обрабоки и выдача сообщений
// */
// struct TrackerGNNAlgorithm {

//     const double CostGate = 9.0;

//     using TypeDetection = Detection<M>;   
//     using TypeEstimator = IMM<M,InitImmFilter1<M>>; 
//     using TypeTrack     = Track<M,TypeEstimator>;
//     using TypeAlgorithm = JVC; 

//     void Bind(TypeDetection & msg) { //const
//         bool isNewTrack = true;

//         //эта операция должна происходить быстро
//         std::vector<TypeTrack*> tracksBindRude;
//         for (auto& track : storageValidTrack) {
//             if (binderRude(track, msg)) { // проверяю привязывается отметка с большим стробом. !!!Какой здесь строб?
//                 tracksBindRude.push_back(&track);
//             }
//         }

//         //тут возможно нужно запустить расчёт парралельно, посмотри std::future!!!
//         for (auto& track : tracksBindRude) {
//             double distance = binder(*track, msg); // считаем дистанцию точно
//             if (distance < CostGate) { // что за costgate?
//                 association(*track, msg, distance);
//                 isNewTrack = false;
//             }
//         }

//         if (isNewTrack) {
//             auto track = storageTrack.GetTrack();
//             track.Initialization(msg);
//             storageValidTrack.InsertTrack(track);
//             // извлекаю трассу из хранилища инициализирую
//         }
//     }

//     void Step(SectorLabel const& msg) {
//         std::vector<TypeTrack*> trackChange;
//         for (auto& track : storageValidTrack) {
//             if ((association.checkTrackAssociation(track)==false) && (cover(track, msg))) { // если не было привязки, если секторная метка прошла по трассе
//                         //association(track)????
//                 track.step(msg.Time); //<- делаю пустой шаг
//                 trackChange.push_back(&track);
//             }
//         }

//         association.Step(msg, trackChange);//<- выполняю шаг ассоциации
//         former.ToUser(trackChange);//<-выдаю на потребителей если трасса подтвежденная

//         std::vector<TypeTrack*> trackDeleted;
//         for (auto& track : storageValidTrack) { // ищу все трассы без измерений
//             // if (IsDeleted(&track)) { //!!!!!!!!!!!!!!!! ПРОВЕРКА
//                 trackDeleted.push_back(&track);
//             // }
//         }
//         storageValidTrack.Delete(trackDeleted);
//         storageTrack.Insert(trackDeleted);
//         former.ToUserDeleted(trackDeleted);
//     }

// private:
//     BinderRude<TypeDetection, TypeTrack> binderRude;
//     Binder<TypeDetection, TypeTrack> binder;

//     Cover<TypeTrack> cover; //класс проверяющий что секторная метка прошла трассу
//     FormerToUser<TypeTrack> former;// класс выдающий данный на потребителя

//     Association<TypeAlgorithm, TypeDetection, TypeTrack> association;

//     StorageValidTrack<TypeTrack> storageValidTrack; // хранилище с иницилизированными трассами
//     StorageTrack<TypeTrack> storageTrack; // хранилище со всеми трассами

// };

// /**
// * Базовый класс, запуск в потоке и обработка сообщений
// */
// struct TrackerGNN {
//     void Run() {
//         try {
//             for (;;) {
//                 incoming.Wait()
//                     .Handle<Detection<M>> ( 
//                         [this](Detection<M> & msg) { //const
//                             // проходишь по всех трассам и пробуешь привязать
//                             algorithm.Bind(msg);
//                         })
//                     .Handle<SectorLabel> (
//                         [this](SectorLabel const& msg) {
//                             //в обрабобтчике секторных меток нужно какой кластер полностью осмотрен
//                             //и где нужно сделать шаг
//                             algorithm.Step(msg);
//                         });
//                 // добавляешь новые сообщения для обработки если необходимо
//             }
//         } catch(CloseTask const&) {
//             std::cout << "[TrackerGNN] Done " << std::endl;
//         }
//     }

//     void Done() {
//         GetSender().Send(CloseTask());
//     }

//     Channel::Sender GetSender() {
//         return incoming;
//     }

//     Channel::Sender ToCdgParser;

// private:
//     TrackerGNNAlgorithm algorithm;


//     Channel::Receiver incoming;
// };

// static bool terminate_signal = false;

// void term_handler(int) {
//     terminate_signal = true;
// }







// TEST_CASE("proto")
// {

// // int main(int argc, char* argv[]) {
// //     //обработка argc и argv, передача кофигов и прочего

// //     //
//     CdgParser cdgParser;
//     TrackerGNN tracker;

//     //
//     cdgParser.ToTrackerGNN = std::move(tracker.GetSender());
//     tracker.ToCdgParser = std::move(cdgParser.GetSender());

//     //
//     std::thread threadPhyCdg([&cdgParser] {cdgParser.Run();}); // запускаем обработку кодограмм
//     std::thread threadTracking([&tracker] {tracker.Run();}); // математика обработки
//     // и т.д. запускаем все процессы

//     while (!terminate_signal) {
//         auto next_begin
//             = std::chrono::system_clock::now() + std::chrono::milliseconds(10);
//         auto end = std::chrono::system_clock::now();
//         // тут отправляешь KeepAlive и если необходимо отправляешь системной вермя в tracker

//         std::this_thread::sleep_for(
//             std::chrono::duration<float, std::milli>(next_begin - end)
//         );
//     }

//     //
//     tracker.Done();
//     cdgParser.Done();

//     //
//     threadPhyCdg.join();
//     threadTracking.join();
// // }

//     BENCHMARK("proto_bench"){

//     };
// }