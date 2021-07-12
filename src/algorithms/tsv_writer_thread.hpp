
#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <thread>
#include "atomic_queue.h"
#include <limits>

typedef std::numeric_limits< double > dbl;

namespace odgi {
	namespace algorithms {

		class tsv_writer {

		private:
			// BED files are 0-based http://genome.ucsc.edu/FAQ/FAQformat#format1
			struct tsv_record_t {
				const std::string query;
				const std::string target;
				bool walk_from_front;
			};

			std::thread writer_thread;
			atomic_queue::AtomicQueue2<tsv_record_t*, 2 << 16> tsv_record_queue;
			std::atomic<bool> work_todo;
			ofstream out_stream;

		public:

			tsv_writer(void) {
				work_todo.store(false);
			}

			~tsv_writer(void) {
				close_writer();
			}

			void writer_func(void) {
				auto* tsv_record = new tsv_record_t();
				while (work_todo.load() || !tsv_record_queue.was_empty()) {
					if (tsv_record_queue.try_pop(tsv_record)) {
						do {
							// BED files are 0-based http://genome.ucsc.edu/FAQ/FAQformat#format1
							out_stream.precision(dbl::max_digits10);
							out_stream << tsv_record->query << "\t" << tsv_record->target << "\t" << tsv_record->walk_from_front << std::endl;
						} while (tsv_record_queue.try_pop(tsv_record));
					} else {
						std::this_thread::sleep_for(std::chrono::nanoseconds(1));
					}
				}
			}

			// open backing file and start writer_thread
			void open_writer(std::string o_file) {
				if (!work_todo.load()) {
					out_stream = ofstream(o_file);
					work_todo.store(true);
					writer_thread = std::thread(&tsv_writer::writer_func, this);
				}
			}

			void close_writer(void) {
				if (work_todo.load()) {
					work_todo.store(false);
					while (!tsv_record_queue.was_empty()) {
						std::this_thread::sleep_for(std::chrono::milliseconds(1));
					}
					if (writer_thread.joinable()) {
						writer_thread.join();
					}
				}
				if (out_stream.is_open()) {
					out_stream.close();
				}
			}

			/// write into our write buffer
			/// open_writer() must be called first to set up our buffer and writer
			void append(const std::string& query, const std::string& target, const bool walk_from_front) {
				tsv_record_queue.push(new tsv_record_t{query, target, walk_from_front});
			}
		};

	}

}