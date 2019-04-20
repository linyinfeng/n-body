#ifndef N_BODY_OVERLOADED_HPP
#define N_BODY_OVERLOADED_HPP

namespace n_body::overloaded {

template <class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template <class... Ts> overloaded(Ts...)->overloaded<Ts...>;

} // namespace n_body::overloaded

#endif
