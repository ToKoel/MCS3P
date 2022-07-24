//
//  ProgressBar.hpp
//  MCS3P
//
//  Created by Tobias Köhler on 13.05.21.
//

#ifndef ProgressBar_h
#define ProgressBar_h

#include <atomic>
#include <mutex>
#include <iostream>
#include <string>

class ProgressBar{
public:
    ProgressBar(){initialize();};
    
    void step(double value, const std::string& status){
        set_progress(value);
        set_status_text(status);
    }
    
    void initialize(){
        set_bar_width(50);
        fill_bar_progress_with("■");
        fill_bar_remainder_with(" ");
        update(0);
    }
    
    void set_progress(double value){
        //std::unique_lock<std::mutex> lock{mutex_};
        progress_ = value;
    }
    
    void set_bar_width(size_t width) {
        std::unique_lock<std::mutex> lock{mutex_};
        bar_width_ = width;
    }

    void fill_bar_progress_with(const std::string& chars) {
        std::unique_lock<std::mutex> lock{mutex_};
        fill_ = chars;
    }

    void fill_bar_remainder_with(const std::string& chars) {
        std::unique_lock<std::mutex> lock{mutex_};
        remainder_ = chars;
    }

    void set_status_text(const std::string& status) {
        std::unique_lock<std::mutex> lock{mutex_};
        status_text_ = status;
    }
    
    void update(float value, std::ostream &os = std::cout) {
        set_progress(value);
        write_progress(os);
    }

    void write_progress(std::ostream &os = std::cout) {
        std::unique_lock<std::mutex> lock{mutex_};

        // No need to write once progress is 100%
        if (progress_ > 100.0f) return;

        // move up one line
        os << "\033[A";
        // erase current line
        os << "\33[2K";
        // Move cursor to the first position on the same line and flush
        os << "\r" << std::flush;

        // Start bar
        os << "[";

        const auto completed = static_cast<size_t>(progress_ * static_cast<float>(bar_width_));
        for (size_t i = 0; i < bar_width_; ++i) {
          if (i <= completed)
            os << fill_;
          else
            os << remainder_;
        }

        // End bar
        os << "]";

        // Write progress percentage
        os << " " << std::min(static_cast<size_t>(progress_*100.0), size_t(100)) << "%";

        // Write status text
        os << " " << status_text_;
        
        os << "\n";
    }
    
private:
    std::mutex mutex_;
    double progress_{0.0f};
    size_t bar_width_{60};
    std::string fill_{"#"}, remainder_{" "}, status_text_{""};
};


#endif /* ProgressBar_h */
