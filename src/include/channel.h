#pragma once

#include <mutex>
#include <condition_variable>
#include <queue>
#include <memory>

namespace Channel {

struct IMessage {
    virtual ~IMessage() {}
};

template <class TypeValue>
struct WrapperMessage : public IMessage {
    TypeValue Value;
    WrapperMessage (TypeValue const& val)
        : Value(val) {}
};

struct Queue {

    template <class TypeMsg>
    void Push(TypeMsg const& msg) {
        std::lock_guard<std::mutex> lg(m);
        q.push(std::make_shared<WrapperMessage<TypeMsg>>(msg));
        c.notify_all();
    }

    void Push(std::shared_ptr<IMessage>& msg) {
        std::lock_guard<std::mutex> lg(m);
        q.push(msg);
        c.notify_all();
    }

    std::shared_ptr<IMessage> WaitAndPop() {
        std::unique_lock<std::mutex> ul(m);
        c.wait(ul, [this] {
            return !this->q.empty();
        });
        auto res = q.front();
        q.pop();
        return res;
    }

private:
    std::mutex m;
    std::condition_variable c;
    std::queue<std::shared_ptr<IMessage>> q;
};

struct Sender {
    //Sender() = delete;
    Sender& operator= (Sender const&) = delete;

    Sender(Sender&& other) {
        if (&other != this) {
            q = other.q;
            other.q = nullptr;
        }
    }

    Sender& operator= (Sender&& other) {
        if (&other == this) {
            return *this;
        }
        q = other.q;
        other.q = nullptr;
        return *this;
    }


    explicit Sender()
    {}

    explicit Sender(Queue *q_) : q(q_)
    {}

    template <typename Message>
    void Send (const Message& msg) {
        if (q != nullptr) {
            q->Push(msg);
        }
    }

    void Send(std::shared_ptr<IMessage>& msg) {
        if (q != nullptr) {
            q->Push(msg);
        }
    }

private:
    Queue *q = nullptr;
};

template <typename PreviousDispatcher,
          typename Msg,
          typename Func>
class TemplateDispatcher {

    Queue* q;
    PreviousDispatcher *prev;
    Func f;
    bool chained;

    TemplateDispatcher(TemplateDispatcher const&) = delete;
    TemplateDispatcher& operator= (TemplateDispatcher const&) = delete;


    template <typename Dispatcher,
              typename OtherMsg,
              typename OtherFunc>
    friend class TemplateDispatcher;

    void waitAndDispatch() {
        for (;;) {
            auto msg = q->WaitAndPop();
            if (dispatch(msg)) {
                break;
            }
        }
    }

    bool dispatch(std::shared_ptr<IMessage> const& msg) {
        if (WrapperMessage<Msg> *wrapper=
                dynamic_cast<WrapperMessage<Msg>*>(msg.get())) {
            f(wrapper->Value);
            return true;
        } else {
            return prev->dispatch(msg);
        }
    }

public:

    TemplateDispatcher(TemplateDispatcher&& other):
        q(other.q),
        prev(other.prev),
        f(std::move(f)),
        chained(other.chained) {
        other.chained = true;
    }

    TemplateDispatcher(Queue* q_,
                       PreviousDispatcher* prev_,
                       Func&& f_):
        q(q_), prev(prev_), f(std::forward<Func>(f_)),
        chained(false)
    {
        prev_->chained = true;
    }


    template <typename OtherMsg,
              typename OtherFunc>
    TemplateDispatcher<TemplateDispatcher,
                       OtherMsg,
                       OtherFunc> Handle(OtherFunc&& of) {
        return TemplateDispatcher<TemplateDispatcher,
                                  OtherMsg,
                                  OtherFunc>(q, this, std::forward<OtherFunc>(of));
    }
    ~TemplateDispatcher() noexcept(false) {
        if (!chained) {
            waitAndDispatch();
        }
    }
};


struct CloseTask
{};

struct Dispatcher {
private:
    Queue* q;
    bool chained;

    Dispatcher(Dispatcher const&) = delete;
    Dispatcher& operator= (Dispatcher const&) = delete;

    template <typename Dispatcher,
              typename Msg,
              typename Func>
    friend struct TemplateDispatcher;
    void waitAndDispatch() {
        for (;;) {
            auto msg = q->WaitAndPop();
            dispatch(msg);
        }
    }
    bool dispatch (std::shared_ptr<IMessage> const& msg) {
        if (dynamic_cast<WrapperMessage<CloseTask>*>(msg.get())) {
            throw CloseTask();
        }
        return false;
    }

public:

    Dispatcher(Dispatcher&& other):
        q(other.q),
        chained(other.chained) {
        other.chained = true;
    }
    explicit Dispatcher(Queue* q_) :
        q(q_), chained(false) {}

    template <typename Message, typename Func>
    TemplateDispatcher<Dispatcher, Message, Func>
    Handle (Func&& f) {
        return TemplateDispatcher<Dispatcher, Message, Func>(
                    q, this, std::forward<Func>(f));
    }

    ~Dispatcher() noexcept(false) {
        if (!chained) {
            waitAndDispatch();
        }
    }
};

struct Receiver {
    operator Sender() {
        return Sender(&q);
    }
    Dispatcher Wait() {
        return Dispatcher(&q);
    }
private:
    Queue q;
};

}
